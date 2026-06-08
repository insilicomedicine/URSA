from __future__ import annotations

import argparse
import gzip
import json
import logging
import sys
from pathlib import Path

from retrocast import adapt_routes
from retrocast import chem
from retrocast import get_adapter
from retrocast import Target
from retrocast.exceptions import ChemError
from retrocast.io import load_collected_routes

from ursa.basic.node import RetrosyntheticNode
from ursa.basic.path import RetrosyntheticPath
from ursa.errors import InvalidTargetKeyError

_BUILTIN_BENCHMARKS = ("EXPERT_2026", "DRUGS_CLINICALS_2026", "USPTO_190")


class _CliFormatter(logging.Formatter):
    """Human-readable formatter for ursa logger _payload records."""

    def format(self, record: logging.LogRecord) -> str:
        """Render ``record._payload`` as a one-line human-readable CLI message.

        Supports three payload shapes: per-path summaries (``is_route_solved``
        present), dataset summaries (``solv_2`` present), and arbitrary
        fallback payloads rendered via ``str``.

        :param record: Log record produced by
            :class:`~ursa.logging.Logger`.
        :type record: logging.LogRecord

        :return: A single formatted line suitable for terminal output.
        :rtype: str
        """
        payload: dict = getattr(record, "_payload", {})
        if "is_route_solved" in payload:
            mark = "+" if payload["is_route_solved"] else "-"
            score = payload.get("chemcensor_per_route", 0.0)
            steps = payload.get("num_steps", "?")
            pid = payload.get("path_id", "")
            consistent = "" if payload.get("is_consistent", True) else " [inconsistent]"
            bb = "" if payload.get("all_bb_found", True) else " [bb missing]"
            return f"  [{mark}] {pid}  score={score:.2f}  steps={steps}{consistent}{bb}"
        if "solv_2" in payload:
            m = payload
            return (
                f"  Solv-2={m['solv_2']:.4f}  "
                f"solved={m['solved_routes']}/{m['total_molecules']}  "
                f"mean_score={m['mean_chemcensor_score']:.4f}  "
                f"passed_steps={m['passed_steps']}/{m['total_steps']}"
            )
        return str(payload)


def _setup_logging() -> None:
    """Attach a stderr handler with :class:`_CliFormatter` to the ursa logger.

    Configures ``ursa.ursa`` at ``INFO`` level so that pipeline events
    produced by :class:`~ursa.logging.Logger` are printed in a compact,
    human-readable form during a CLI run.
    """
    handler = logging.StreamHandler(sys.stderr)
    handler.setFormatter(_CliFormatter())
    ursa_logger = logging.getLogger("ursa.ursa")
    ursa_logger.addHandler(handler)
    ursa_logger.setLevel(logging.INFO)


def _build_parser() -> argparse.ArgumentParser:
    """Build the ``ursa-bench`` argument parser.

    Defines the mutually exclusive input sources (``--routes`` vs
    ``--input``/``--adapter``), benchmark selection, output directory,
    and optional CSV column names for custom benchmarks.

    :return: Configured :class:`argparse.ArgumentParser`.
    :rtype: argparse.ArgumentParser
    """
    p = argparse.ArgumentParser(
        prog="ursa-bench",
        description="Score retrosynthesis predictions with URSA.",
    )

    src = p.add_mutually_exclusive_group(required=True)
    src.add_argument(
        "--routes",
        metavar="ROUTES_JSON_GZ",
        help="Already-adapted RetroCast routes (.json.gz).",
    )
    src.add_argument(
        "--input",
        metavar="RAW_JSON_GZ",
        help="Raw model predictions (.json.gz) — requires --adapter.",
    )

    p.add_argument(
        "--adapter",
        metavar="ADAPTER",
        help=(
            "Adapter name for --input (e.g. retrochimera, aizynth, dms). "
            "Run ``retrocast list-adapters`` to see all options."
        ),
    )
    p.add_argument(
        "--benchmark",
        metavar="BENCHMARK",
        required=True,
        help=(
            f"Built-in preset ({', '.join(_BUILTIN_BENCHMARKS)}) "
            "or path to a CSV file. "
            "When a CSV path is given, --id-col and --smiles-col can be used "
            "to specify column names (defaults: 'Structure ID', 'SMILES')."
        ),
    )
    p.add_argument(
        "--output",
        metavar="OUTPUT_DIR",
        required=True,
        help="Directory to write metrics.json and best_paths.json.",
    )
    p.add_argument(
        "--stem",
        metavar="STEM",
        default="",
        help="Optional filename prefix (e.g. 'retrochimera_expert').",
    )
    p.add_argument(
        "--id-col",
        default="Structure ID",
        metavar="COL",
        help="ID column name for custom benchmark CSV (default: 'Structure ID').",
    )
    p.add_argument(
        "--smiles-col",
        default="SMILES",
        metavar="COL",
        help="SMILES column name for custom benchmark CSV (default: 'SMILES').",
    )
    return p


def _load_routes(routes_path: str) -> dict:
    """Load already-adapted RetroCast routes from a ``.json.gz`` file.

    Expects a *collected routes* archive (schema 2) as written by
    :func:`retrocast.io.save_collected_routes`: a JSON object mapping each
    target id to a list of :class:`~retrocast.Route` objects.

    :param routes_path: Path to a routes archive produced by RetroCast.
    :type routes_path: str

    :return: RetroCast routes dict keyed by target identifier.
    :rtype: dict
    """
    return load_collected_routes(Path(routes_path))


def _make_target(key: str) -> Target | None:
    """Build a RetroCast :class:`~retrocast.Target` from a raw input key.

    The raw predictions file is keyed by target SMILES, so the key is used
    as the target id while its canonical form is used as the SMILES. The
    InChIKey is derived from the SMILES. Returns ``None`` when the key is
    not a valid SMILES, in which case the adapter falls back to deriving
    the target from each raw route.

    :param key: Target identifier / SMILES from the raw predictions file.
    :type key: str

    :return: A populated :class:`~retrocast.Target`, or ``None``.
    :rtype: retrocast.Target | None
    """
    try:
        try:
            return Target(
                id=key,
                smiles=chem.canonicalize_smiles(key),
                inchikey=chem.get_inchi_key(key),
            )
        except ChemError as exc:
            raise InvalidTargetKeyError(key) from exc
    except InvalidTargetKeyError:
        return None


def _adapt_and_load(raw_path: str, adapter_name: str) -> dict:
    """Adapt raw model predictions into RetroCast routes and group them.

    Reads a gzipped JSON object keyed by target SMILES/ID and runs the
    named RetroCast adapter over each entry via
    :func:`retrocast.adapt_routes`, which validates and canonicalises the
    routes. Returns the adapted routes grouped by their original key.

    :param raw_path: Path to the raw ``.json.gz`` predictions file.
    :type raw_path: str
    :param adapter_name: Name of the RetroCast adapter to apply
        (e.g. ``"retrochimera"``, ``"aizynth"``, ``"dms"``).
    :type adapter_name: str

    :return: RetroCast routes dict ready for :func:`_routes_to_paths`.
    :rtype: dict
    """

    adapter = get_adapter(adapter_name)

    with gzip.open(raw_path, "rb") as f:
        raw_data = json.loads(f.read())

    if not isinstance(raw_data, dict):
        print(
            "error: --input file must be a JSON object keyed by target SMILES or ID",
            file=sys.stderr,
        )
        sys.exit(1)

    routes: dict = {}
    for key, payload in raw_data.items():
        adapted = adapt_routes(payload, adapter, target=_make_target(str(key)))
        if adapted:
            routes[key] = adapted

    return routes


def _routes_to_paths(routes_dict: dict):
    """Convert a RetroCast routes dict into URSA :class:`RetrosyntheticPath` objects.

    Each target produces one path per ranked route; path ids follow the
    ``"{target_id}__r{rank}"`` convention.

    :param routes_dict: RetroCast routes keyed by target id.
    :type routes_dict: dict

    :return: List of :class:`~ursa.RetrosyntheticPath` instances.
    :rtype: list[RetrosyntheticPath]
    """

    def mol_to_node(mol):
        """Recursively convert a RetroCast molecule into a :class:`RetrosyntheticNode`.

        :param mol: RetroCast molecule node with an optional
            ``product_of`` reaction containing ``reactants``.
        :return: Fully populated :class:`RetrosyntheticNode`.
        :rtype: RetrosyntheticNode
        """
        if mol.product_of is None:
            return RetrosyntheticNode(smiles=mol.smiles)
        children = tuple(mol_to_node(r) for r in mol.product_of.reactants)
        return RetrosyntheticNode(smiles=mol.smiles, children=children)

    paths = []
    for target_id, route_list in routes_dict.items():
        for rank, route in enumerate(route_list):
            path_id = f"{target_id}__r{rank + 1}"
            root = mol_to_node(route.target)
            paths.append(RetrosyntheticPath(path_id=path_id, root=root))
    return paths


def _load_benchmark(benchmark: str, id_col: str, smiles_col: str):
    """Resolve a benchmark name or CSV path to a :class:`BenchmarkDataset`.

    Built-in names (``EXPERT_2026``, ``DRUGS_CLINICALS_2026``,
    ``USPTO_190``) are returned as the bundled presets; anything else is
    treated as a CSV path and loaded via
    :meth:`BenchmarkDataset.from_csv`. Exits with code 1 if the value is
    neither a known preset nor an existing file.

    :param benchmark: Preset name or filesystem path to a CSV file.
    :type benchmark: str
    :param id_col: Column name for molecule identifiers (CSV path only).
    :type id_col: str
    :param smiles_col: Column name for SMILES strings (CSV path only).
    :type smiles_col: str

    :return: Loaded :class:`BenchmarkDataset`.
    :rtype: BenchmarkDataset
    """
    from ursa import BenchmarkDataset

    if benchmark == "EXPERT_2026":
        return BenchmarkDataset.EXPERT_2026
    if benchmark == "DRUGS_CLINICALS_2026":
        return BenchmarkDataset.DRUGS_CLINICALS_2026
    if benchmark == "USPTO_190":
        return BenchmarkDataset.USPTO_190

    csv_path = Path(benchmark)
    if not csv_path.exists():
        print(
            f"error: benchmark not recognised and file not found: {benchmark}",
            file=sys.stderr,
        )
        sys.exit(1)
    return BenchmarkDataset.from_csv(csv_path, id_col=id_col, smiles_col=smiles_col)


def main(argv: list[str] | None = None) -> None:
    """Entry point for ``ursa-bench``."""
    parser = _build_parser()
    args = parser.parse_args(argv)

    if args.input and not args.adapter:
        parser.error("--adapter is required when using --input")

    _setup_logging()

    print("Loading routes...")
    if args.routes:
        routes_dict = _load_routes(args.routes)
    else:
        routes_dict = _adapt_and_load(args.input, args.adapter)

    paths = _routes_to_paths(routes_dict)
    print(f"  {len(paths)} paths loaded for {len(routes_dict)} targets")

    dataset = _load_benchmark(args.benchmark, args.id_col, args.smiles_col)
    target_smiles = list(dataset.target_smiles)
    print(f"Benchmark: {args.benchmark} ({len(target_smiles)} targets)")

    from ursa import Ursa

    print("Scoring...")
    ursa = Ursa()
    result = ursa.score_dataset(paths, target_smiles=target_smiles)

    m = result.metrics
    print(
        f"\nResults:\n"
        f"  total_molecules:       {m.total_molecules}\n"
        f"  molecules_with_route:  {m.molecules_with_route}\n"
        f"  solved_routes:         {m.solved_routes}\n"
        f"  Solv-2:                {m.solv_2:.4f}\n"
        f"  mean_chemcensor_score: {m.mean_chemcensor_score:.4f}\n"
        f"  passed_steps:          {m.passed_steps} / {m.total_steps}"
    )

    metrics_path, paths_path = result.save(args.output, stem=args.stem)
    print(f"\nSaved:\n  {metrics_path}\n  {paths_path}")


if __name__ == "__main__":
    main()
