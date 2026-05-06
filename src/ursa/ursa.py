from __future__ import annotations

from os import PathLike
from typing import Sequence

from rdkit import Chem

from .basic import DatasetResult
from .basic import PathResult
from .basic import RetrosyntheticPath
from .collapsing import PathCollapser
from .configs import DataConfig
from .logging import get_logger
from .logging import Logger
from .metrics import DatasetMetricsCalculator
from .scoring import BestPathSelector
from .scoring import PathScorer
from .scoring import ReactionScorer
from .validation import BuildingBlockChecker
from .validation import PathConsistencyChecker


class Ursa:
    """Public API for retrosynthetic path evaluation.

    Orchestrates the full evaluation pipeline:

    1. **Consistency check** — verifies tree structure via
       :class:`~ursa.PathConsistencyChecker`.
    2. **Building-block check** — looks up leaf molecules in the catalog
       via :class:`~ursa.BuildingBlockChecker`.
    3. **Collapsing** — generates all valid collapsed variants via
       :class:`~ursa.PathCollapser`.
    4. **Scoring** — scores each variant via :class:`~ursa.PathScorer`.
    5. **Selection** — picks the best variant via
       :class:`~ursa.BestPathSelector`.
    6. **Metrics** — aggregates dataset-level metrics via
       :class:`~ursa.DatasetMetricsCalculator` (``score_dataset`` only).

    The concrete reaction scorer is supplied by the caller, keeping
    :class:`Ursa` decoupled from any specific scoring implementation.

    :param scorer: Reaction scorer implementing
        :class:`~ursa.ReactionScorer`. The canonical choice is
        ``chemcensor.ChemCensor``. When ``None``, a
        ``chemcensor.ChemCensor`` is instantiated with
        :attr:`~ursa.DataConfig.chemcensor_db_path`.
    :type scorer: ReactionScorer | None
    :param bb_catalog_path: Path to the building-block catalog file
        (plain-text SMILES or CSV with a ``smiles`` column). Defaults to
        :attr:`~ursa.DataConfig.bb_catalog_path`.
    :type bb_catalog_path: str | PathLike | None
    """

    def __init__(
        self,
        *,
        scorer: ReactionScorer | None = None,
        bb_catalog_path: str | PathLike | None = None,
    ) -> None:
        """Initialize Ursa.

        Constructs all internal pipeline components. Both parameters are
        optional and fall back to the bundled defaults in
        :class:`~ursa.DataConfig`.

        :param scorer: Reaction scorer implementing
            :class:`~ursa.ReactionScorer`. When ``None``, a
            ``chemcensor.ChemCensor`` is created from the default DB.
        :type scorer: ReactionScorer | None
        :param bb_catalog_path: Path to the building-block catalog.
            When ``None``, uses :attr:`~ursa.DataConfig.bb_catalog_path`.
        :type bb_catalog_path: str | PathLike | None
        """
        if scorer is None:
            from chemcensor import ChemCensor

            scorer = ChemCensor(db_path=DataConfig.chemcensor_db_path)
        if bb_catalog_path is None:
            bb_catalog_path = DataConfig.bb_catalog_path
        self._consistency_checker = PathConsistencyChecker()
        self._bb_checker = BuildingBlockChecker.from_file(bb_catalog_path)
        self._collapser = PathCollapser()
        self._path_scorer = PathScorer(scorer)
        self._selector = BestPathSelector()
        self._metrics_calculator = DatasetMetricsCalculator()
        self._logger: Logger = get_logger(__name__)

    def score(self, path: RetrosyntheticPath) -> PathResult:
        """Evaluate a single retrosynthetic path.

        Runs the full pipeline (consistency → BB check → collapse →
        score → select) and returns a :class:`~ursa.PathResult` for the
        best-scoring collapsed variant.

        :param path: The retrosynthetic path to evaluate.
        :type path: RetrosyntheticPath

        :return: Full evaluation result including consistency, building
            blocks, and best-variant scoring.
        :rtype: PathResult
        """
        is_consistent = self._consistency_checker.check(path)
        starting_materials = self._bb_checker.check(path)
        all_bb_found = all(bb.found_in_catalog for bb in starting_materials)

        variants = self._collapser.collapse(path)
        variant_results = tuple(self._path_scorer.score(v) for v in variants)
        best_variant = self._selector.select(variant_results)

        is_route_solved = (
            is_consistent and all_bb_found and best_variant.all_steps_passed
        )

        result = PathResult(
            original_path=path,
            is_consistent=is_consistent,
            starting_materials=starting_materials,
            all_bb_found=all_bb_found,
            best_variant=best_variant,
            is_route_solved=is_route_solved,
        )
        self._logger.log_path(result)
        return result

    def score_dataset(
        self,
        paths: Sequence[RetrosyntheticPath],
        target_smiles: Sequence[str] | None = None,
        total_molecules: int | None = None,
    ) -> DatasetResult:
        """Evaluate a dataset of retrosynthetic paths.

        When ``target_smiles`` is provided, only paths whose root molecule
        matches a target are evaluated. Multiple paths for the same target
        are scored and the best :class:`~ursa.PathResult` is kept (solved
        routes preferred, ties broken by variant quality). Targets with no
        matching path count toward ``total_molecules`` but produce no
        :class:`~ursa.PathResult`, so they lower ``solv_2``.

        When ``target_smiles`` is omitted, every path is scored as-is and
        ``total_molecules`` defaults to ``len(paths)``.

        :param paths: Retrosynthetic paths to evaluate.
        :type paths: Sequence[RetrosyntheticPath]
        :param target_smiles: Canonical input-target SMILES. Paths whose
            root is not in this set are discarded. Mutually exclusive with
            ``total_molecules``.
        :type target_smiles: Sequence[str] | None
        :param total_molecules: Denominator for ``solv_2``.
            Ignored when ``target_smiles`` is given. Defaults to
            ``len(paths)`` when neither argument is supplied.
        :type total_molecules: int | None

        :raises ValueError: If both ``target_smiles`` and
            ``total_molecules`` are provided.

        :return: Per-path results and aggregate dataset metrics.
        :rtype: DatasetResult
        """
        if target_smiles is not None and total_molecules is not None:
            raise ValueError(
                "Provide either target_smiles or total_molecules, not both."
            )

        if target_smiles is not None:
            groups = self._group_paths_by_target(paths, target_smiles)
            total_molecules = len(groups)
            scored: list[PathResult] = []
            for target_paths in groups.values():
                if not target_paths:
                    continue
                results = tuple(self.score(p) for p in target_paths)
                best = (
                    self._select_best_path_result(results)
                    if len(results) > 1
                    else results[0]
                )
                scored.append(best)
            path_results = tuple(scored)
        else:
            if total_molecules is None:
                total_molecules = len(paths)
            path_results = tuple(self.score(path) for path in paths)

        metrics = self._metrics_calculator.calculate(path_results, total_molecules)
        dataset_result = DatasetResult(path_results=path_results, metrics=metrics)
        self._logger.log_dataset(dataset_result)
        return dataset_result

    def _group_paths_by_target(
        self,
        paths: Sequence[RetrosyntheticPath],
        target_smiles: Sequence[str],
    ) -> dict[str, list[RetrosyntheticPath]]:
        """Group paths by canonical target SMILES, discarding off-target paths.

        Preserves the order of unique canonical targets as they appear in
        ``target_smiles``. Invalid or duplicate target SMILES are silently
        skipped.

        :param paths: All candidate paths to filter and group.
        :param target_smiles: Requested target SMILES strings.
        :return: Ordered mapping from canonical SMILES to matching paths.
        """
        canonical_targets: dict[str, None] = {}
        for smi in target_smiles:
            mol = Chem.MolFromSmiles(smi)
            if mol is not None:
                canonical_targets[Chem.MolToSmiles(mol)] = None

        groups: dict[str, list[RetrosyntheticPath]] = {
            cs: [] for cs in canonical_targets
        }
        for path in paths:
            cs = path.root.canonical_smiles
            if cs in groups:
                groups[cs].append(path)
        return groups

    def _select_best_path_result(
        self,
        path_results: Sequence[PathResult],
    ) -> PathResult:
        """Select the best :class:`~ursa.PathResult` from multiple candidates.

        Solved routes are preferred over unsolved. Within each group,
        the candidate whose ``best_variant`` ranks highest under
        :class:`~ursa.BestPathSelector` is returned.

        :param path_results: Non-empty sequence of scored path results for
            the same target molecule.
        :return: Best :class:`~ursa.PathResult` among the candidates.
        """
        solved = [pr for pr in path_results if pr.is_route_solved]
        pool = solved if solved else list(path_results)
        best_vr = self._selector.select(tuple(pr.best_variant for pr in pool))
        return next(pr for pr in pool if pr.best_variant is best_vr)
