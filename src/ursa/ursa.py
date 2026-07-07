from __future__ import annotations

import logging
from dataclasses import dataclass
from os import PathLike
from typing import Sequence

from rdkit import Chem

from .basic import BuildingBlock
from .basic import DatasetResult
from .basic import PathResult
from .basic import RetrosyntheticPath
from .basic import VariantResult
from .collapsing import PathCollapser
from .configs import DataConfig
from .logging import get_logger
from .logging import Logger
from .metrics import DatasetMetricsCalculator
from .scoring import BestPathSelector
from .scoring import PathScorer
from .scoring import ReactionScore
from .scoring import ReactionScorer
from .validation import BuildingBlockChecker
from .validation import PathConsistencyChecker

# Plain (message-only) progress logger; the CLI handler renders these as
# human-readable lines. Distinct from the structured per-event ``Logger``.
_progress = logging.getLogger(__name__)

# Emit a prepare/assemble progress line at most every this many paths.
_PROGRESS_EVERY = 50


@dataclass(frozen=True)
class _PreparedPath:
    """Scorer-independent checks for a path, computed *before* collapsing.

    Solv-0 (stock termination) depends only on these original-path
    properties, never on collapsed variants. Paths that fail Solv-0
    cannot pass Solv-1/2 either, so collapsing and reaction scoring are
    skipped for them entirely (see :meth:`Ursa._score_paths_parallel`).
    """

    path: RetrosyntheticPath
    is_consistent: bool
    starting_materials: tuple[BuildingBlock, ...]
    all_bb_found: bool
    is_no_synthesis: bool

    @property
    def passes_solv_0(self) -> bool:
        """``True`` if the original path satisfies stock termination (Solv-0)."""
        return not self.is_no_synthesis and self.is_consistent and self.all_bb_found


class Ursa:
    """Public API for retrosynthetic path evaluation.

    Orchestrates the full evaluation pipeline:

    1. **Consistency check** — verifies tree structure via
       :class:`~ursa.PathConsistencyChecker`.
    2. **Building-block check** — looks up leaf molecules in the catalog
       via :class:`~ursa.BuildingBlockChecker`. Together with step 1 this
       decides Solv-0; paths failing Solv-0 skip steps 3-5 entirely.
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
        parallel: bool = False,
        n_workers: int | None = None,
        chemcensor_db_path: str | PathLike | None = None,
    ) -> None:
        """Initialize Ursa.

        Constructs all internal pipeline components. All parameters are
        optional and fall back to the bundled defaults in
        :class:`~ursa.DataConfig`.

        :param scorer: Reaction scorer implementing
            :class:`~ursa.ReactionScorer`. When ``None``, a
            ``chemcensor.ChemCensor`` is created from the default DB.
        :type scorer: ReactionScorer | None
        :param bb_catalog_path: Path to the building-block catalog.
            When ``None``, uses :attr:`~ursa.DataConfig.bb_catalog_path`.
        :type bb_catalog_path: str | PathLike | None
        :param parallel: When ``True``, :meth:`score_dataset` scores all
            reactions in parallel via
            :func:`chemcensor.parallel.score_batch` (multiprocess). This
            bypasses ``scorer`` and uses ``chemcensor_db_path`` directly.
            :meth:`score` (single path) is always sequential.
        :type parallel: bool
        :param n_workers: Number of scorer worker processes for parallel
            mode. ``None`` auto-scales to the host CPU count.
        :type n_workers: int | None
        :param chemcensor_db_path: Explicit path to a ChemCensor
            ``.sqlite`` database. Used by both the default sequential
            scorer and parallel mode. When ``None``, the default database
            is downloaded on first use from Hugging Face (see
            :func:`~ursa.data.chemcensor_db.ensure_chemcensor_db`) to
            :attr:`~ursa.DataConfig.chemcensor_db_path`.
        :type chemcensor_db_path: str | PathLike | None
        """
        if chemcensor_db_path is None and (scorer is None or parallel):
            from .data.chemcensor_db import ensure_chemcensor_db

            chemcensor_db_path = ensure_chemcensor_db(DataConfig.chemcensor_db_path)
        self._chemcensor_db_path = chemcensor_db_path
        if scorer is None:
            from chemcensor import ChemCensor

            scorer = ChemCensor(db_path=self._chemcensor_db_path)
        if bb_catalog_path is None:
            from .data.bb_catalog import ensure_bb_catalog

            bb_catalog_path = ensure_bb_catalog(DataConfig.bb_catalog_path)
        self._consistency_checker = PathConsistencyChecker()
        self._bb_checker = BuildingBlockChecker.from_file(bb_catalog_path)
        self._collapser = PathCollapser()
        self._path_scorer = PathScorer(scorer)
        self._selector = BestPathSelector()
        self._metrics_calculator = DatasetMetricsCalculator()
        self._logger: Logger = get_logger(__name__)
        self._parallel = parallel
        self._n_workers = n_workers

    def score(self, path: RetrosyntheticPath) -> PathResult:
        """Evaluate a single retrosynthetic path.

        Runs the full pipeline (consistency → BB check → collapse →
        score) and selects the best collapsed variant **independently for
        each Solv level**. Solv-0 (stock termination) is variant
        independent: it holds when the tree is consistent and every
        starting material is in the catalog. Solv-1 / Solv-2 additionally
        require *some* variant whose every step passes ChemCensor without
        / with functional-group matching.

        :param path: The retrosynthetic path to evaluate.
        :type path: RetrosyntheticPath

        :return: Full evaluation result with per-level pass flags and the
            best variant selected for each Solv level.
        :rtype: PathResult
        """
        prepared = self._prepare(path)
        if not prepared.passes_solv_0:
            return self._finalize_unscored(prepared)
        variants = self._collapser.collapse(path)
        variant_results = self._path_scorer.score_variants(variants)
        return self._finalize(prepared, variant_results)

    def _prepare(self, path: RetrosyntheticPath) -> _PreparedPath:
        """Run the cheap, scorer-independent checks for ``path``.

        Performs only the consistency check and building-block lookup —
        everything Solv-0 depends on. Collapsing and reaction scoring are
        deliberately *not* done here: they are expensive and only matter
        for paths that already pass Solv-0.

        :param path: The retrosynthetic path to prepare.
        :type path: RetrosyntheticPath

        :return: Bundle of Solv-0 pre-scoring results.
        :rtype: _PreparedPath
        """
        is_consistent = self._consistency_checker.check(path)
        starting_materials = self._bb_checker.check(path)
        all_bb_found = all(bb.found_in_catalog for bb in starting_materials)
        return _PreparedPath(
            path=path,
            is_consistent=is_consistent,
            starting_materials=starting_materials,
            all_bb_found=all_bb_found,
            is_no_synthesis=path.num_steps == 0,
        )

    def _finalize(
        self,
        prepared: _PreparedPath,
        variant_results: tuple[VariantResult, ...],
    ) -> PathResult:
        """Select per-level best variants and assemble the path result.

        :param prepared: Output of :meth:`_prepare` for the path.
        :type prepared: _PreparedPath
        :param variant_results: Scored variants (original + collapsed) of
            the prepared path.
        :type variant_results: tuple[VariantResult, ...]

        :return: Full evaluation result for the path (also logged).
        :rtype: PathResult
        """
        best_solv_1 = self._selector.select_for_solv_1(variant_results)
        best_solv_2 = self._selector.select_for_solv_2(variant_results)

        passes_solv_0 = prepared.passes_solv_0
        passes_solv_1 = passes_solv_0 and best_solv_1.all_steps_pass_solv_1
        passes_solv_2 = passes_solv_0 and best_solv_2.all_steps_pass_solv_2

        result = PathResult(
            original_path=prepared.path,
            is_consistent=prepared.is_consistent,
            starting_materials=prepared.starting_materials,
            all_bb_found=prepared.all_bb_found,
            is_no_synthesis=prepared.is_no_synthesis,
            passes_solv_0=passes_solv_0,
            passes_solv_1=passes_solv_1,
            passes_solv_2=passes_solv_2,
            best_variant_solv_1=best_solv_1,
            best_variant_solv_2=best_solv_2,
        )
        self._logger.log_path(result)
        return result

    def _finalize_unscored(self, prepared: _PreparedPath) -> PathResult:
        """Assemble the result for a path that fails Solv-0.

        Such a path cannot pass any Solv level, so its reactions are never
        collapsed or scored. The best variants are empty placeholders over
        the original path (no step results), which keeps them out of the
        dataset mean-score diagnostics.

        :param prepared: Output of :meth:`_prepare` for the path.
        :type prepared: _PreparedPath

        :return: Full evaluation result for the path (also logged).
        :rtype: PathResult
        """
        empty = VariantResult(
            path=prepared.path,
            step_results=(),
            all_steps_pass_solv_1=False,
            all_steps_pass_solv_2=False,
            mean_score_without_fg=0.0,
            mean_score_with_fg=0.0,
        )
        result = PathResult(
            original_path=prepared.path,
            is_consistent=prepared.is_consistent,
            starting_materials=prepared.starting_materials,
            all_bb_found=prepared.all_bb_found,
            is_no_synthesis=prepared.is_no_synthesis,
            passes_solv_0=False,
            passes_solv_1=False,
            passes_solv_2=False,
            best_variant_solv_1=empty,
            best_variant_solv_2=empty,
        )
        self._logger.log_path(result)
        return result

    def _score_many(
        self, paths: Sequence[RetrosyntheticPath]
    ) -> tuple[PathResult, ...]:
        """Score a sequence of paths, parallelising when enabled.

        :param paths: Paths to evaluate.
        :type paths: Sequence[RetrosyntheticPath]

        :return: One :class:`~ursa.PathResult` per input path, in order.
        :rtype: tuple[PathResult, ...]
        """
        if self._parallel and paths:
            return self._score_paths_parallel(paths)
        return tuple(self.score(path) for path in paths)

    def _score_paths_parallel(
        self, paths: Sequence[RetrosyntheticPath]
    ) -> tuple[PathResult, ...]:
        """Score many paths with a single bulk, multiprocess scoring pass.

        Runs the cheap Solv-0 checks for every path, then collapses and
        scores **only** the paths that pass Solv-0 (paths failing Solv-0
        cannot pass Solv-1/2, so their reactions are never scored). The
        unique reactions across the surviving variants are scored once via
        :func:`chemcensor.parallel.score_batch`, then each path result is
        assembled from the shared score map.

        :param paths: Paths to evaluate.
        :type paths: Sequence[RetrosyntheticPath]

        :return: One :class:`~ursa.PathResult` per input path, in order.
        :rtype: tuple[PathResult, ...]
        """
        total = len(paths)
        _progress.info("Preparing %d paths (consistency, building blocks)...", total)
        prepared = [self._prepare(path) for path in paths]

        scored_idx = [i for i, item in enumerate(prepared) if item.passes_solv_0]
        _progress.info(
            "  %d/%d paths pass Solv-0; collapsing those...", len(scored_idx), total
        )
        variants_by_idx: dict[int, tuple[RetrosyntheticPath, ...]] = {}
        for count, i in enumerate(scored_idx, start=1):
            variants_by_idx[i] = self._collapser.collapse(prepared[i].path)
            if count % _PROGRESS_EVERY == 0 or count == len(scored_idx):
                _progress.info("  collapsed %d/%d paths", count, len(scored_idx))

        unique_reactions: dict[str, None] = {}
        for variants in variants_by_idx.values():
            for variant in variants:
                for node in variant.get_all_steps():
                    unique_reactions.setdefault(node.reaction_smiles, None)

        _progress.info(
            "Scoring %d unique reactions (%s workers)...",
            len(unique_reactions),
            self._n_workers or "auto",
        )
        raw_scores = self._parallel_score_reactions(list(unique_reactions))

        _progress.info("Assembling %d path results...", total)
        results = []
        for i, item in enumerate(prepared):
            if i in variants_by_idx:
                variant_results = self._path_scorer.build_results_from_raw(
                    variants_by_idx[i], raw_scores
                )
                results.append(self._finalize(item, variant_results))
            else:
                results.append(self._finalize_unscored(item))
        return tuple(results)

    def _parallel_score_reactions(
        self, reactions: Sequence[str]
    ) -> dict[str, ReactionScore]:
        """Score unique reactions in parallel via ChemCensor's pipeline.

        :param reactions: Distinct forward reaction SMILES to score.
        :type reactions: Sequence[str]

        :return: Mapping from reaction SMILES to its (raw) score result.
        :rtype: dict[str, ReactionScore]
        """
        if not reactions:
            return {}
        from chemcensor.parallel import ParallelConfig
        from chemcensor.parallel import score_batch

        config = ParallelConfig(n_workers=self._n_workers)
        results = score_batch(
            reactions, db_path=self._chemcensor_db_path, config=config
        )
        return dict(zip(reactions, results))

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
            flat_paths = [p for target_paths in groups.values() for p in target_paths]
            flat_results = self._score_many(flat_paths)
            result_by_path = {id(p): r for p, r in zip(flat_paths, flat_results)}
            scored: list[PathResult] = []
            for target_paths in groups.values():
                if not target_paths:
                    continue
                results = tuple(result_by_path[id(p)] for p in target_paths)
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
            path_results = self._score_many(list(paths))

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

        Candidates are ranked by the highest Solv level they satisfy
        (Solv-2 > Solv-1 > Solv-0), breaking ties by the mean
        functional-group-aware score of the Solv-2 best variant.

        :param path_results: Non-empty sequence of scored path results for
            the same target molecule.
        :return: Best :class:`~ursa.PathResult` among the candidates.
        """
        return max(
            path_results,
            key=lambda pr: (
                pr.passes_solv_2,
                pr.passes_solv_1,
                pr.passes_solv_0,
                pr.best_variant_solv_2.mean_score_with_fg,
            ),
        )
