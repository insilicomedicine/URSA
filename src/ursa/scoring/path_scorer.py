from typing import Sequence

from ..basic import RetrosyntheticNode
from ..basic import RetrosyntheticPath
from ..basic import StepResult
from ..basic import VariantResult
from ..configs import PathScoringConfig
from .base import ReactionScorer


class PathScorer:
    """Scores retrosynthetic steps for one or many path variants.
    For a single variant, use :meth:`score`. For multiple collapsed variants,
    use :meth:`score_variants` to score each unique ``reaction_smiles`` once
    and build results from a shared cache.
    """

    def __init__(self, scorer: ReactionScorer) -> None:
        """Initialize PathScorer.

        :param scorer: Reaction scorer implementing
            :class:`~ursa.ReactionScorer`.
        :type scorer: ReactionScorer
        """
        self._scorer = scorer

    def score(self, path: RetrosyntheticPath) -> VariantResult:
        """Score all internal nodes of ``path`` and return a
        :class:`~ursa.VariantResult`.

        Traverses the tree in DFS pre-order. For each internal node,
        looks up its reaction score. Aggregates ``percent_found``,
        ``chemcensor_per_route``, and ``all_steps_passed`` from the
        individual step results.

        :param path: The retrosynthetic path variant to score.
        :type path: RetrosyntheticPath

        :return: Aggregated scoring result for the variant.
        :rtype: VariantResult
        """
        scores = self._score_unique_reactions((path,))
        return self._build_variant_result(path, scores)

    def score_variants(
        self, variants: Sequence[RetrosyntheticPath]
    ) -> tuple[VariantResult, ...]:
        """Score multiple path variants, scoring each unique reaction once.

        Collapsed variants of the same path (produced by
        :class:`~ursa.PathCollapser`) share many identical reactions.
        This method first collects the set of unique ``reaction_smiles``
        across **all** ``variants``, scores each exactly once, and then
        builds every :class:`~ursa.VariantResult` from that shared cache.
        This avoids the redundant re-scoring of shared reactions that
        would occur when scoring each variant in isolation, which is the
        dominant cost for deep, many-step paths.

        :param variants: Path variants to score (e.g. the output of
            :meth:`~ursa.PathCollapser.collapse`).
        :type variants: Sequence[RetrosyntheticPath]

        :return: One :class:`~ursa.VariantResult` per input variant, in
            the same order.
        :rtype: tuple[VariantResult, ...]
        """
        scores = self._score_unique_reactions(variants)
        return tuple(self._build_variant_result(v, scores) for v in variants)

    def _score_unique_reactions(
        self, variants: Sequence[RetrosyntheticPath]
    ) -> dict[str, float]:
        """Collect every distinct reaction across ``variants`` and score each once.

        The returned mapping holds the *normalised* score for each
        reaction (negative raw scores are clamped to
        :attr:`PathScoringConfig.failed_step_score`). Reactions are
        scored in first-seen DFS order so call order is deterministic.

        :param variants: Path variants whose internal nodes are scored.
        :type variants: Sequence[RetrosyntheticPath]

        :return: Mapping from ``reaction_smiles`` to its normalised score.
        :rtype: dict[str, float]
        """
        unique_reactions: dict[str, None] = {}
        for variant in variants:
            for node in variant.get_all_steps():
                unique_reactions.setdefault(node.reaction_smiles, None)

        return {
            reaction: self._normalise_score(self._scorer.score(reaction))
            for reaction in unique_reactions
        }

    def _build_variant_result(
        self, path: RetrosyntheticPath, scores: dict[str, float]
    ) -> VariantResult:
        """Assemble a :class:`~ursa.VariantResult` for ``path`` from cached scores.

        :param path: The path variant to build a result for.
        :type path: RetrosyntheticPath
        :param scores: Mapping from ``reaction_smiles`` to normalised
            score, as produced by :meth:`_score_unique_reactions`. Must
            contain every reaction reachable from ``path``.
        :type scores: dict[str, float]

        :return: Aggregated scoring result for the variant.
        :rtype: VariantResult
        """
        step_results = tuple(
            self._build_step_result(n, scores) for n in path.get_all_steps()
        )
        total = len(step_results)

        if total == 0:
            return VariantResult(
                path=path,
                step_results=(),
                percent_found=0.0,
                chemcensor_per_route=0.0,
                all_steps_passed=True,
            )

        passed_count = sum(1 for sr in step_results if sr.passed)
        percent_found = passed_count / total
        chemcensor_per_route = sum(sr.score for sr in step_results) / total
        all_steps_passed = all(sr.passed for sr in step_results)

        return VariantResult(
            path=path,
            step_results=step_results,
            percent_found=percent_found,
            chemcensor_per_route=chemcensor_per_route,
            all_steps_passed=all_steps_passed,
        )

    def _build_step_result(
        self, node: RetrosyntheticNode, scores: dict[str, float]
    ) -> StepResult:
        """Construct a :class:`~ursa.StepResult` for ``node`` from cached scores.

        A step is considered *passed* when its score is at or above
        :attr:`PathScoringConfig.pass_threshold`.

        :param node: Internal node whose ``reaction_smiles`` score is
            looked up.
        :type node: RetrosyntheticNode
        :param scores: Mapping from ``reaction_smiles`` to normalised
            score; must contain ``node.reaction_smiles``.
        :type scores: dict[str, float]

        :return: Scoring result for this step.
        :rtype: StepResult
        """
        score = scores[node.reaction_smiles]
        passed = score >= PathScoringConfig.pass_threshold.value
        return StepResult(node=node, score=score, passed=passed)

    @staticmethod
    def _normalise_score(raw: float) -> float:
        """Clamp a raw scorer value, mapping negative (failure) values to a floor.

        :param raw: Raw score returned by the reaction scorer. Negative
            values indicate scoring failures.
        :type raw: float

        :return: ``PathScoringConfig.failed_step_score`` if ``raw`` is
            negative, otherwise ``float(raw)``.
        :rtype: float
        """
        return PathScoringConfig.failed_step_score.value if raw < 0 else float(raw)
