from __future__ import annotations

import json
from dataclasses import dataclass
from pathlib import Path
from typing import Any

from ..configs import PathScoringConfig
from .building_block import BuildingBlock
from .node import RetrosyntheticNode
from .path import RetrosyntheticPath


@dataclass(frozen=True)
class StepResult:
    """Scoring result for a single retrosynthetic step.

    Holds both ChemCensor score variants for the step's reaction. The
    functional-group-agnostic score drives Solv-1, the
    functional-group-aware score drives Solv-2. A step *passes* a level
    when the corresponding score is strictly greater than
    :attr:`PathScoringConfig.pass_threshold`.

    :param node: The internal tree node whose reaction was scored.
    :type node: RetrosyntheticNode
    :param score_without_fg: Score ignoring functional-group matching
        (reaction-center presence only). Negative scorer failures are
        normalised to ``0.0`` before being stored here.
    :type score_without_fg: float
    :param score_with_fg: Score requiring functional-group matching.
        Normalised the same way as ``score_without_fg``.
    :type score_with_fg: float
    """

    node: RetrosyntheticNode
    score_without_fg: float
    score_with_fg: float

    @property
    def passes_solv_1(self) -> bool:
        """``True`` if the FG-agnostic score exceeds ``pass_threshold``."""
        return self.score_without_fg > PathScoringConfig.pass_threshold.value

    @property
    def passes_solv_2(self) -> bool:
        """``True`` if the FG-aware score exceeds ``pass_threshold``."""
        return self.score_with_fg > PathScoringConfig.pass_threshold.value

    def to_dict(self) -> dict[str, Any]:
        """Serialise to a plain dict with both scores and pass flags."""
        return {
            "smiles": self.node.smiles,
            "reaction_smiles": self.node.reaction_smiles,
            "score_without_fg": self.score_without_fg,
            "score_with_fg": self.score_with_fg,
            "passes_solv_1": self.passes_solv_1,
            "passes_solv_2": self.passes_solv_2,
        }


@dataclass(frozen=True)
class VariantResult:
    """Scoring result for one collapsed variant of a retrosynthetic path.

    A variant is obtained by :class:`~ursa.PathCollapser` and differs
    from the original path by having some non-adjacent steps merged.

    :param path: The (possibly collapsed) retrosynthetic path variant.
    :type path: RetrosyntheticPath
    :param step_results: Per-step scoring results in DFS pre-order.
    :type step_results: tuple[StepResult, ...]
    :param all_steps_pass_solv_1: ``True`` if the variant has at least one
        step and every step passes Solv-1 (``score_without_fg > 0``).
    :type all_steps_pass_solv_1: bool
    :param all_steps_pass_solv_2: ``True`` if the variant has at least one
        step and every step passes Solv-2 (``score_with_fg > 0``).
    :type all_steps_pass_solv_2: bool
    :param mean_score_without_fg: Mean functional-group-agnostic score
        across all steps (``0.0`` for an empty variant).
    :type mean_score_without_fg: float
    :param mean_score_with_fg: Mean functional-group-aware score across
        all steps (``0.0`` for an empty variant).
    :type mean_score_with_fg: float
    """

    path: RetrosyntheticPath
    step_results: tuple[StepResult, ...]
    all_steps_pass_solv_1: bool
    all_steps_pass_solv_2: bool
    mean_score_without_fg: float
    mean_score_with_fg: float

    def to_dict(self) -> dict[str, Any]:
        """Serialise to a plain dict with per-step scores and aggregates."""
        return {
            "num_steps": self.path.num_steps,
            "all_steps_pass_solv_1": self.all_steps_pass_solv_1,
            "all_steps_pass_solv_2": self.all_steps_pass_solv_2,
            "mean_score_without_fg": self.mean_score_without_fg,
            "mean_score_with_fg": self.mean_score_with_fg,
            "steps": [sr.to_dict() for sr in self.step_results],
        }


@dataclass(frozen=True)
class PathResult:
    """Full evaluation result for a single retrosynthetic path.

    Combines consistency, building-block availability, and the best
    collapsed variant **selected independently for each Solv level**.

    The Solv hierarchy is:

    * **Solv-0** (stock termination): valid tree (``is_consistent``) and
      all starting materials in the catalog (``all_bb_found``).
    * **Solv-1**: Solv-0 plus some variant where every step passes
      ChemCensor *without* functional-group matching.
    * **Solv-2**: Solv-0 plus some variant where every step passes
      ChemCensor *with* functional-group matching.

    A path whose target is itself a building block (no reaction steps) is
    flagged ``is_no_synthesis`` and does not pass any Solv level.

    :param original_path: The original (uncollapsed) retrosynthetic path.
    :type original_path: RetrosyntheticPath
    :param is_consistent: ``True`` if the tree passes the consistency check.
    :type is_consistent: bool
    :param starting_materials: Catalog-check results for each leaf node.
    :type starting_materials: tuple[BuildingBlock, ...]
    :param all_bb_found: ``True`` if every starting material is in the catalog.
    :type all_bb_found: bool
    :param is_no_synthesis: ``True`` if the path has no reaction steps.
    :type is_no_synthesis: bool
    :param passes_solv_0: ``True`` if the route satisfies Solv-0.
    :type passes_solv_0: bool
    :param passes_solv_1: ``True`` if the route satisfies Solv-1.
    :type passes_solv_1: bool
    :param passes_solv_2: ``True`` if the route satisfies Solv-2.
    :type passes_solv_2: bool
    :param best_variant_solv_1: Variant minimising Solv-1 failed steps
        (the representative used to decide ``passes_solv_1``).
    :type best_variant_solv_1: VariantResult
    :param best_variant_solv_2: Variant minimising Solv-2 failed steps;
        also the variant used for display.
    :type best_variant_solv_2: VariantResult
    """

    original_path: RetrosyntheticPath
    is_consistent: bool
    starting_materials: tuple[BuildingBlock, ...]
    all_bb_found: bool
    is_no_synthesis: bool
    passes_solv_0: bool
    passes_solv_1: bool
    passes_solv_2: bool
    best_variant_solv_1: VariantResult
    best_variant_solv_2: VariantResult

    def to_dict(self) -> dict[str, Any]:
        """Serialise to a plain dict including both per-level best variants."""
        return {
            "path_id": self.original_path.path_id,
            "is_consistent": self.is_consistent,
            "all_bb_found": self.all_bb_found,
            "is_no_synthesis": self.is_no_synthesis,
            "passes_solv_0": self.passes_solv_0,
            "passes_solv_1": self.passes_solv_1,
            "passes_solv_2": self.passes_solv_2,
            "best_variant_solv_1": self.best_variant_solv_1.to_dict(),
            "best_variant_solv_2": self.best_variant_solv_2.to_dict(),
        }


@dataclass(frozen=True)
class DatasetMetrics:
    """Aggregate metrics computed over a dataset of retrosynthetic paths.

    :param total_molecules: Total number of target molecules in the input
        dataset (including those without a route). Denominator of every
        ``solv_N`` rate.
    :type total_molecules: int
    :param molecules_with_route: Number of targets for which at least one
        retrosynthetic path was provided.
    :type molecules_with_route: int
    :param routes_no_synthesis: Number of evaluated paths with no reaction
        steps (target is itself a building block).
    :type routes_no_synthesis: int
    :param routes_solv_0: Number of routes satisfying Solv-0.
    :type routes_solv_0: int
    :param routes_solv_1: Number of routes satisfying Solv-1.
    :type routes_solv_1: int
    :param routes_solv_2: Number of routes satisfying Solv-2.
    :type routes_solv_2: int
    :param solv_0: ``routes_solv_0 / total_molecules``.
    :type solv_0: float
    :param solv_1: ``routes_solv_1 / total_molecules``.
    :type solv_1: float
    :param solv_2: ``routes_solv_2 / total_molecules``.
    :type solv_2: float
    :param mean_score_without_fg: Mean functional-group-agnostic step
        score across the Solv-1 best variant of every route. Routes that
        fail Solv-0 are not scored (empty best variants) and so do not
        contribute to this mean.
    :type mean_score_without_fg: float
    :param mean_score_with_fg: Mean functional-group-aware step score
        across the Solv-2 best variant of every route. Routes that fail
        Solv-0 do not contribute to this mean.
    :type mean_score_with_fg: float
    """

    total_molecules: int
    molecules_with_route: int
    routes_no_synthesis: int
    routes_solv_0: int
    routes_solv_1: int
    routes_solv_2: int
    solv_0: float
    solv_1: float
    solv_2: float
    mean_score_without_fg: float
    mean_score_with_fg: float


@dataclass(frozen=True)
class DatasetResult:
    """Evaluation result for a full dataset of retrosynthetic paths.

    :param path_results: Per-path evaluation results, one entry per
        evaluated :class:`RetrosyntheticPath`.
    :type path_results: tuple[PathResult, ...]
    :param metrics: Aggregate metrics computed over all ``path_results``.
    :type metrics: DatasetMetrics
    """

    path_results: tuple[PathResult, ...]
    metrics: DatasetMetrics

    def to_dict(self) -> dict[str, Any]:
        """Serialise to a plain dict with ``metrics`` and ``path_results``."""
        m = self.metrics
        return {
            "metrics": {
                "total_molecules": m.total_molecules,
                "molecules_with_route": m.molecules_with_route,
                "routes_no_synthesis": m.routes_no_synthesis,
                "routes_solv_0": m.routes_solv_0,
                "routes_solv_1": m.routes_solv_1,
                "routes_solv_2": m.routes_solv_2,
                "solv_0": m.solv_0,
                "solv_1": m.solv_1,
                "solv_2": m.solv_2,
                "mean_score_without_fg": m.mean_score_without_fg,
                "mean_score_with_fg": m.mean_score_with_fg,
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
