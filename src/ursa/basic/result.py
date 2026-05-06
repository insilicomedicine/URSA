from __future__ import annotations

import json
from dataclasses import dataclass
from pathlib import Path
from typing import Any

from .building_block import BuildingBlock
from .node import RetrosyntheticNode
from .path import RetrosyntheticPath


@dataclass(frozen=True)
class StepResult:
    """Scoring result for a single retrosynthetic step.

    :param node: The tree node (internal) whose reaction was scored.
    :type node: RetrosyntheticNode
    :param score: Raw score returned by the reaction scorer.
    :type score: float
    :param passed: ``True`` if ``score`` is at or above the configured
        pass threshold.
    :type passed: bool
    """

    node: RetrosyntheticNode
    score: float
    passed: bool

    def to_dict(self) -> dict[str, Any]:
        """Serialise to a plain dict with ``smiles``, ``reaction_smiles``,
        ``score``, and ``passed``."""
        return {
            "smiles": self.node.smiles,
            "reaction_smiles": self.node.reaction_smiles,
            "score": self.score,
            "passed": self.passed,
        }


@dataclass(frozen=True)
class VariantResult:
    """Scoring result for one collapsed variant of a retrosynthetic path.

    A variant is obtained by :class:`~ursa.PathCollapser` and differs
    from the original path by having some consecutive steps merged.

    :param path: The (possibly collapsed) retrosynthetic path variant.
    :type path: RetrosyntheticPath
    :param step_results: Scoring results for each internal node of
        ``path``, in DFS pre-order.
    :type step_results: tuple[StepResult, ...]
    :param percent_found: Fraction of steps with ``score > 0``.
    :type percent_found: float
    :param chemcensor_per_route: Mean scorer value across all steps in
        this variant.
    :type chemcensor_per_route: float
    :param all_steps_passed: ``True`` if every step in ``step_results``
        has ``passed == True``.
    :type all_steps_passed: bool
    """

    path: RetrosyntheticPath
    step_results: tuple[StepResult, ...]
    percent_found: float
    chemcensor_per_route: float
    all_steps_passed: bool


@dataclass(frozen=True)
class PathResult:
    """Full evaluation result for a single retrosynthetic path.

    Combines consistency check, building-block availability, and the
    best-scoring collapsed variant selected by
    :class:`~ursa.BestPathSelector`.

    :param original_path: The original (uncollapsed) retrosynthetic path.
    :type original_path: RetrosyntheticPath
    :param is_consistent: ``True`` if the path tree passes the
        consistency check (products match reactants across all steps).
    :type is_consistent: bool
    :param starting_materials: Catalog-check results for each leaf
        (starting material) node of ``original_path``, produced by
        :class:`~ursa.BuildingBlockChecker`.
    :type starting_materials: tuple[BuildingBlock, ...]
    :param all_bb_found: ``True`` if every starting material was found
        in the building-block catalog.
    :type all_bb_found: bool
    :param best_variant: The variant selected by
        :class:`~ursa.BestPathSelector` as the best-scoring collapsed
        form of the path.
    :type best_variant: VariantResult
    :param is_route_solved: ``True`` if the route is fully solved:
        ``is_consistent``, ``all_bb_found``, and
        ``best_variant.all_steps_passed`` are all ``True``.
    :type is_route_solved: bool
    """

    original_path: RetrosyntheticPath
    is_consistent: bool
    starting_materials: tuple[BuildingBlock, ...]
    all_bb_found: bool
    best_variant: VariantResult
    is_route_solved: bool

    def to_dict(self) -> dict[str, Any]:
        """Serialise to a plain dict including ``best_variant`` with per-step
        scores and ``path_id``, ``is_consistent``, ``all_bb_found``,
        ``is_route_solved``."""
        bv = self.best_variant
        return {
            "path_id": self.original_path.path_id,
            "is_consistent": self.is_consistent,
            "all_bb_found": self.all_bb_found,
            "is_route_solved": self.is_route_solved,
            "best_variant": {
                "num_steps": bv.path.num_steps,
                "percent_found": bv.percent_found,
                "chemcensor_per_route": bv.chemcensor_per_route,
                "all_steps_passed": bv.all_steps_passed,
                "steps": [sr.to_dict() for sr in bv.step_results],
            },
        }


@dataclass(frozen=True)
class DatasetMetrics:
    """Aggregate metrics computed over a dataset of retrosynthetic paths.

    :param total_molecules: Total number of molecules in the input
        dataset (including those without a route).
    :type total_molecules: int
    :param molecules_with_route: Number of molecules for which at least
        one retrosynthetic path was provided.
    :type molecules_with_route: int
    :param solved_routes: Number of routes with ``is_route_solved == True``.
    :type solved_routes: int
    :param passed_steps: Total number of steps with ``score > 0`` across
        all best variants.
    :type passed_steps: int
    :param total_steps: Total number of steps across all best variants.
    :type total_steps: int
    :param solv_2: Fraction of solved routes over the total number of
        input molecules: ``solved_routes / total_molecules``.
    :type solv_2: float
    :param mean_chemcensor_score: Mean scorer value across all steps in
        all best variants: ``sum(all scores) / total_steps``.
    :type mean_chemcensor_score: float
    """

    total_molecules: int
    molecules_with_route: int
    solved_routes: int
    passed_steps: int
    total_steps: int
    solv_2: float
    mean_chemcensor_score: float


@dataclass(frozen=True)
class DatasetResult:
    """Evaluation result for a full dataset of retrosynthetic paths.

    :param path_results: Per-path evaluation results, one entry per
        input :class:`RetrosyntheticPath`.
    :type path_results: tuple[PathResult, ...]
    :param metrics: Aggregate metrics computed over all ``path_results``.
    :type metrics: DatasetMetrics
    """

    path_results: tuple[PathResult, ...]
    metrics: DatasetMetrics

    def to_dict(self) -> dict[str, Any]:
        """Serialise to a plain dict with ``metrics`` and ``path_results``
        (each entry is a :meth:`PathResult.to_dict` dump)."""
        m = self.metrics
        return {
            "metrics": {
                "total_molecules": m.total_molecules,
                "molecules_with_route": m.molecules_with_route,
                "solved_routes": m.solved_routes,
                "passed_steps": m.passed_steps,
                "total_steps": m.total_steps,
                "solv_2": m.solv_2,
                "mean_chemcensor_score": m.mean_chemcensor_score,
            },
            "path_results": [pr.to_dict() for pr in self.path_results],
        }

    def save(self, output_dir: Path | str, stem: str = "") -> tuple[Path, Path]:
        """Write metrics and best paths to two JSON files in ``output_dir``.

        Creates ``output_dir`` if it does not exist. File names are
        ``{stem}_metrics.json`` / ``{stem}_best_paths.json`` when *stem* is
        given, otherwise ``metrics.json`` / ``best_paths.json``.

        :param output_dir: Directory to write files into.
        :type output_dir: Path | str
        :param stem: Optional filename prefix.
        :type stem: str

        :return: Tuple of (metrics_path, best_paths_path).
        :rtype: tuple[Path, Path]
        """
        out = Path(output_dir)
        out.mkdir(parents=True, exist_ok=True)
        prefix = f"{stem}_" if stem else ""
        d = self.to_dict()
        metrics_path = out / f"{prefix}metrics.json"
        paths_path = out / f"{prefix}best_paths.json"
        metrics_path.write_text(json.dumps(d["metrics"], indent=2))
        paths_path.write_text(json.dumps(d["path_results"], indent=2))
        return metrics_path, paths_path
