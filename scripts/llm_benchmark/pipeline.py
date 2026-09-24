from __future__ import annotations

import argparse
import csv
import json
import os
import random
import sys
from pathlib import Path

from .sampling import sample_completions
from .scoring import score_records
from .templates import build_prompt

_PRESETS = ("EXPERT_2026", "DRUGS_CLINICALS_2026")
_DEFAULT_ENV_FILE = Path(__file__).resolve().parents[2] / ".env"


def canonicalize_targets(targets: list[str]) -> list[str]:
    """Canonicalize target SMILES with RDKit while preserving stereochemistry."""
    from rdkit import Chem

    canonical = []
    for index, smiles in enumerate(targets, 1):
        molecule = Chem.MolFromSmiles(smiles)
        if molecule is None:
            raise SystemExit(
                f"error: invalid target SMILES at benchmark row {index}: {smiles}"
            )
        canonical.append(
            Chem.MolToSmiles(molecule, canonical=True, isomericSmiles=True)
        )
    return canonical


def load_targets(benchmark: str, smiles_col: str) -> list[str]:
    """Load and canonicalize target SMILES from a preset or custom CSV."""
    if benchmark in _PRESETS:
        from ursa import BenchmarkDataset

        dataset = getattr(BenchmarkDataset, benchmark)
        return canonicalize_targets([entry.smiles for entry in dataset.load()])

    path = Path(benchmark)
    if not path.exists():
        raise SystemExit(
            f"error: benchmark not recognised and file not found: {benchmark}"
        )
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle)
        if reader.fieldnames is None or smiles_col not in reader.fieldnames:
            raise SystemExit(
                f"error: column {smiles_col!r} not in {path}. "
                "Pass --smiles-col (URSA preset CSVs use 'SMILES')."
            )
        targets = [
            value for row in reader if (value := (row.get(smiles_col) or "").strip())
        ]
    if not targets:
        raise SystemExit(f"error: no target SMILES in {path}")
    return canonicalize_targets(targets)


def load_records(path: Path) -> list[dict]:
    """Load and validate completion records from JSONL."""
    if not path.exists():
        return []
    records = []
    with path.open() as handle:
        for line_no, line in enumerate(handle, 1):
            text = line.strip()
            if not text:
                continue
            record = json.loads(text)
            completion = record.get("completion")
            meta = record.get("meta") or {}
            product = meta.get("product_smiles")
            if not isinstance(completion, str) or not isinstance(product, str):
                raise SystemExit(
                    f"error: {path}:{line_no} needs string 'completion' and "
                    "'meta.product_smiles'"
                )
            records.append(record)
    return records


def _env_int(name: str) -> int | None:
    raw = os.environ.get(name)
    if raw is None or raw == "":
        return None
    return int(raw)


def load_environment(argv: list[str]) -> None:
    """Load an env file before environment-backed CLI defaults are evaluated."""
    env_parser = argparse.ArgumentParser(add_help=False)
    env_parser.add_argument("--env-file", type=Path, default=_DEFAULT_ENV_FILE)
    env_args, _ = env_parser.parse_known_args(argv)

    try:
        from dotenv import load_dotenv
    except ModuleNotFoundError as exc:
        raise ModuleNotFoundError(
            "Environment-file support requires the llm-benchmark extra. "
            "Install it with `uv sync --extra llm-benchmark`."
        ) from exc
    load_dotenv(env_args.env_file, override=False)


def build_parser() -> argparse.ArgumentParser:
    """Build the end-to-end benchmark CLI."""
    parser = argparse.ArgumentParser(
        prog="scripts/benchmark_llm.py",
        description="Sample an LLM on multi-step retrosynthesis and score with URSA.",
        epilog=(
            "LiteLLM model examples: openai/gpt-5, gemini/gemini-2.5-pro, "
            "anthropic/claude-sonnet-4-5, xai/grok-4. "
            "Use --base-url for a custom or local endpoint."
        ),
    )
    parser.add_argument(
        "--env-file",
        type=Path,
        default=_DEFAULT_ENV_FILE,
        help="Environment file to load (default: repository .env).",
    )
    parser.add_argument(
        "--benchmark",
        required=True,
        help=(
            "Built-in preset (EXPERT_2026, DRUGS_CLINICALS_2026) or a CSV path. "
            "This list is the denominator for every Solv-N rate."
        ),
    )
    parser.add_argument(
        "--smiles-col",
        default="product_smiles",
        help="SMILES column for a custom CSV (default: product_smiles).",
    )
    parser.add_argument(
        "--limit",
        type=int,
        default=0,
        help="Use only the first N benchmark targets. 0 keeps all.",
    )
    parser.add_argument(
        "--output",
        type=Path,
        required=True,
        help="Directory for metrics JSON and best-path JSON.",
    )
    parser.add_argument(
        "--completions",
        type=Path,
        default=None,
        help="JSONL of completions. Default: OUTPUT/completions.jsonl.",
    )
    parser.add_argument(
        "--stem", default="llm", help="Filename prefix for score files."
    )
    parser.add_argument(
        "--model",
        default="",
        help="LiteLLM model identifier for this benchmark run.",
    )
    parser.add_argument(
        "--base-url",
        default=(
            os.environ.get("LLM_BASE_URL")
            or os.environ.get("AZURE_OPENAI_ENDPOINT")
            or None
        ),
        help=(
            "Optional provider endpoint. Reads LLM_BASE_URL or "
            "AZURE_OPENAI_ENDPOINT."
        ),
    )
    parser.add_argument(
        "--api-key",
        default=(
            os.environ.get("LLM_API_KEY")
            or os.environ.get("AZURE_OPENAI_API_KEY")
            or None
        ),
        help=(
            "Provider key. Reads LLM_API_KEY or AZURE_OPENAI_API_KEY; optional "
            "when LiteLLM can read a provider-specific variable."
        ),
    )
    parser.add_argument(
        "--api-version",
        default=(
            os.environ.get("LLM_API_VERSION")
            or os.environ.get("AZURE_OPENAI_API_VERSION")
            or os.environ.get("AZURE_API_VERSION")
            or None
        ),
        help="Optional API version, required by some Azure OpenAI deployments.",
    )
    parser.add_argument("--samples", type=int, default=10, help="Samples per target.")
    parser.add_argument(
        "--request-workers",
        "--api-workers",
        dest="request_workers",
        type=int,
        default=8,
        help="Concurrent inference requests (default: 8).",
    )
    parser.add_argument("--timeout", type=float, default=180.0)
    parser.add_argument(
        "--temperature",
        type=float,
        default=None,
        help="Sampling temperature. Omit to use the provider default.",
    )
    top_k_default = _env_int("URSA_TOP_K")
    parser.add_argument(
        "--top-k",
        type=int,
        default=top_k_default if top_k_default is not None else 10,
        help="Routes kept per target after parsing. 0 keeps all. Default: 10.",
    )
    parser.add_argument(
        "--score-workers",
        type=int,
        default=_env_int("URSA_WORKERS"),
        help="ChemCensor worker processes. Omit for sequential; 0 auto-scales.",
    )
    parser.add_argument(
        "--bb-catalog",
        type=Path,
        default=os.environ.get("URSA_BB_CATALOG") or None,
    )
    parser.add_argument(
        "--chemcensor-db",
        type=Path,
        default=os.environ.get("URSA_CHEMCENSOR_DB") or None,
    )
    parser.add_argument(
        "--fresh",
        action="store_true",
        help="Delete an existing completions file before sampling.",
    )
    parser.add_argument("--sample-only", action="store_true")
    parser.add_argument("--score-only", action="store_true")
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Print one prompt for the first target and exit.",
    )
    parser.add_argument("--seed", type=int, default=None)
    return parser


def main(argv: list[str] | None = None) -> None:
    """Run sampling and scoring, or either phase independently."""
    argv = sys.argv[1:] if argv is None else argv
    load_environment(argv)
    parser = build_parser()
    args = parser.parse_args(argv)
    if args.sample_only and args.score_only:
        parser.error("choose one of --sample-only and --score-only")
    if args.samples < 1:
        parser.error("--samples must be >= 1")
    if args.limit < 0:
        parser.error("--limit must be >= 0")
    if args.request_workers < 1:
        parser.error("--request-workers must be >= 1")
    if args.seed is not None:
        random.seed(args.seed)

    targets = load_targets(args.benchmark, args.smiles_col)
    if args.limit:
        targets = targets[: args.limit]
    print(f"Benchmark: {args.benchmark} ({len(targets)} targets)")
    if args.dry_run:
        print(build_prompt(targets[0]))
        return

    completions_path = args.completions or (args.output / "completions.jsonl")
    if not args.score_only:
        if not args.model:
            parser.error("--model is required when sampling")
        sample_completions(
            targets,
            load_records(completions_path),
            n_samples=args.samples,
            completions_path=completions_path,
            base_url=args.base_url,
            api_key=args.api_key,
            api_version=args.api_version,
            model=args.model,
            timeout=args.timeout,
            temperature=args.temperature,
            request_workers=args.request_workers,
            fresh=args.fresh,
        )
        print(f"Completions: {completions_path}")

    if args.sample_only:
        return
    if not completions_path.exists():
        raise SystemExit(f"error: completions file not found: {completions_path}")

    records = load_records(completions_path)
    print(f"Scoring {len(records)} completions")
    score_records(
        records,
        targets,
        top_k=args.top_k,
        score_workers=args.score_workers,
        bb_catalog=args.bb_catalog,
        chemcensor_db=args.chemcensor_db,
        output=args.output,
        stem=args.stem,
    )
