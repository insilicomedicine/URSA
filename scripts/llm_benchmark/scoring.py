from __future__ import annotations

import json
from collections import defaultdict
from pathlib import Path


def _to_node(molecule):
    """Convert a RetroCast molecule tree into an URSA node tree."""
    from ursa import RetrosyntheticNode

    if molecule.product_of is None:
        return RetrosyntheticNode(smiles=molecule.smiles)
    children = tuple(_to_node(reactant) for reactant in molecule.product_of.reactants)
    return RetrosyntheticNode(smiles=molecule.smiles, children=children)


def completions_to_paths(
    records: list[dict], targets: list[str], top_k: int
) -> tuple[list, int, int]:
    """Parse completions and return paths, adapted count, and attempted count."""
    from retrocast import Target
    from retrocast import adapt_routes
    from retrocast import chem
    from retrocast import get_adapter
    from ursa import RetrosyntheticPath

    wanted = set(targets)
    grouped: dict[str, list[dict]] = defaultdict(list)
    for record in records:
        smiles = record["meta"]["product_smiles"]
        if smiles in wanted:
            grouped[smiles].append(record)

    adapter = get_adapter("ursa")
    paths = []
    n_adapted = 0
    n_completions = 0
    for target, group in grouped.items():
        n_completions += len(group)
        target_obj = Target(
            id=target,
            smiles=chem.canonicalize_smiles(target),
            inchikey=chem.get_inchi_key(target),
        )
        routes = adapt_routes(group, adapter, target=target_obj)
        n_adapted += len(routes)
        if top_k > 0:
            routes = routes[:top_k]
        for rank, route in enumerate(routes, 1):
            paths.append(
                RetrosyntheticPath(
                    path_id=f"{target}__r{rank}",
                    root=_to_node(route.target),
                )
            )
    return paths, n_adapted, n_completions


def score_records(
    records: list[dict],
    targets: list[str],
    *,
    top_k: int,
    score_workers: int | None,
    bb_catalog: Path | None,
    chemcensor_db: Path | None,
    output: Path,
    stem: str,
) -> None:
    """Score records while retaining the full target list as denominator."""
    from ursa import Ursa

    paths, n_adapted, n_completions = completions_to_paths(records, targets, top_k)
    rate = (n_adapted / n_completions) if n_completions else 0.0
    print(
        f"Adapted {n_adapted}/{n_completions} completions ({rate:.4f}); "
        f"{len(paths)} routes kept"
    )

    parallel = score_workers is not None
    n_workers = score_workers if score_workers else None
    ursa = Ursa(
        parallel=parallel,
        n_workers=n_workers,
        bb_catalog_path=bb_catalog,
        chemcensor_db_path=chemcensor_db,
    )
    result = ursa.score_dataset(paths, target_smiles=targets)
    metrics = result.metrics
    print(
        f"\nResults:\n"
        f"  total_molecules:       {metrics.total_molecules}\n"
        f"  molecules_with_route:  {metrics.molecules_with_route}\n"
        f"  Solv-0 (STR):          {metrics.solv_0:.4f}  "
        f"({metrics.routes_solv_0} routes)\n"
        f"  Solv-1:                {metrics.solv_1:.4f}  "
        f"({metrics.routes_solv_1} routes)\n"
        f"  Solv-2:                {metrics.solv_2:.4f}  "
        f"({metrics.routes_solv_2} routes)\n"
        f"  mean_score_without_fg: {metrics.mean_score_without_fg:.4f}\n"
        f"  mean_score_with_fg:    {metrics.mean_score_with_fg:.4f}\n"
        f"  adapt_rate:            {rate:.4f}"
    )
    metrics_path, paths_path = result.save(output, stem=stem)
    metrics_payload = json.loads(metrics_path.read_text())
    metrics_payload.update(
        {
            "completions_total": n_completions,
            "completions_adapted": n_adapted,
            "adapt_rate": rate,
            "routes_kept": len(paths),
        }
    )
    metrics_path.write_text(json.dumps(metrics_payload, indent=2) + "\n")
    print(f"\nSaved:\n  {metrics_path}\n  {paths_path}")
