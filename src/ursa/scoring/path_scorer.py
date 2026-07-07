from typing import Mapping
from typing import Sequence

from ..basic import RetrosyntheticNode
from ..basic import RetrosyntheticPath
from ..basic import StepResult
from ..basic import VariantResult
from ..configs import PathScoringConfig
from .base import ReactionScore
from .base import ReactionScorer


class PathScorer:
    """Scores retrosynthetic steps for one or many path variants.

    For each internal node the injected :class:`~ursa.ReactionScorer` is
    asked to :meth:`~ursa.ReactionScorer.evaluate` the reaction, yielding
    both the functional-group-agnostic (Solv-1) and functional-group-aware
    (Solv-2) scores in a single call. Use :meth:`score` for a single
    variant and :meth:`score_variants` for a batch of collapsed variants
    sharing reactions.

    :param scorer: The reaction scorer used to evaluate individual steps.
    :type scorer: ReactionScorer
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
        """Score multiple path variants, evaluating each unique reaction once.

        Collapsed variants of the same path share many identical
        reactions. This method collects the set of unique
        ``reaction_smiles`` across **all** ``variants``, evaluates each
        exactly once, and builds every :class:`~ursa.VariantResult` from
        that shared cache.

        :param variants: Path variants to score (e.g. the output of
            :meth:`~ursa.PathCollapser.collapse`).
        :type variants: Sequence[RetrosyntheticPath]

        :return: One :class:`~ursa.VariantResult` per input variant, in
            the same order.
        :rtype: tuple[VariantResult, ...]
        """
        scores = self._score_unique_reactions(variants)
        return tuple(self._build_variant_result(v, scores) for v in variants)

    def build_results_from_raw(
        self,
        variants: Sequence[RetrosyntheticPath],
        raw_scores: Mapping[str, ReactionScore],
    ) -> tuple[VariantResult, ...]:
        """Build variant results from an externally computed score map.

        Used when reactions are scored in bulk outside the scorer (e.g.
        by :func:`chemcensor.parallel.score_batch`). ``raw_scores`` must
        cover every ``reaction_smiles`` reachable from ``variants``; its
        values are normalised here exactly as :meth:`score_variants` does.

        :param variants: Path variants to assemble results for.
        :type variants: Sequence[RetrosyntheticPath]
        :param raw_scores: Mapping from ``reaction_smiles`` to a
            :class:`~ursa.ReactionScore` (with raw, un-normalised scores).
        :type raw_scores: Mapping[str, ReactionScore]

        :return: One :class:`~ursa.VariantResult` per input variant.
        :rtype: tuple[VariantResult, ...]
        """
        scores = {
            reaction: (
                self._normalise_score(result.without_functional_groups),
                self._normalise_score(result.with_functional_groups),
            )
            for reaction, result in raw_scores.items()
        }
        return tuple(self._build_variant_result(v, scores) for v in variants)

    def _score_unique_reactions(
        self, variants: Sequence[RetrosyntheticPath]
    ) -> dict[str, tuple[float, float]]:
        """Evaluate every distinct reaction across ``variants`` once.

        The returned mapping holds the *normalised* ``(without_fg,
        with_fg)`` scores for each reaction (negative scorer failures are
        clamped to :attr:`PathScoringConfig.failed_step_score`). Reactions
        are evaluated in first-seen DFS order for deterministic call order.

        :param variants: Path variants whose internal nodes are scored.
        :type variants: Sequence[RetrosyntheticPath]

        :return: Mapping from ``reaction_smiles`` to ``(without_fg, with_fg)``.
        :rtype: dict[str, tuple[float, float]]
        """
        unique_reactions: dict[str, None] = {}
        for variant in variants:
            for node in variant.get_all_steps():
                unique_reactions.setdefault(node.reaction_smiles, None)

        scores: dict[str, tuple[float, float]] = {}
        for reaction in unique_reactions:
            result = self._scorer.evaluate(reaction)
            scores[reaction] = (
                self._normalise_score(result.without_functional_groups),
                self._normalise_score(result.with_functional_groups),
            )
        return scores

    def _build_variant_result(
        self,
        path: RetrosyntheticPath,
        scores: dict[str, tuple[float, float]],
    ) -> VariantResult:
        """Assemble a :class:`~ursa.VariantResult` for ``path`` from cached scores.

        :param path: The path variant to build a result for.
        :type path: RetrosyntheticPath
        :param scores: Mapping from ``reaction_smiles`` to ``(without_fg,
            with_fg)``; must contain every reaction reachable from ``path``.
        :type scores: dict[str, tuple[float, float]]

        :return: Aggregated scoring result for the variant.
        :rtype: VariantResult
        """
        step_results = tuple(
            self._build_step_result(n, scores) for n in path.get_all_steps()
        )
        has_steps = len(step_results) > 0

        if not has_steps:
            return VariantResult(
                path=path,
                step_results=(),
                all_steps_pass_solv_1=False,
                all_steps_pass_solv_2=False,
                mean_score_without_fg=0.0,
                mean_score_with_fg=0.0,
            )

        total = len(step_results)
        return VariantResult(
            path=path,
            step_results=step_results,
            all_steps_pass_solv_1=all(sr.passes_solv_1 for sr in step_results),
            all_steps_pass_solv_2=all(sr.passes_solv_2 for sr in step_results),
            mean_score_without_fg=sum(sr.score_without_fg for sr in step_results)
            / total,
            mean_score_with_fg=sum(sr.score_with_fg for sr in step_results) / total,
        )

    def _build_step_result(
        self,
        node: RetrosyntheticNode,
        scores: dict[str, tuple[float, float]],
    ) -> StepResult:
        """Construct a :class:`~ursa.StepResult` for ``node`` from cached scores.

        :param node: Internal node whose ``reaction_smiles`` score is
            looked up.
        :type node: RetrosyntheticNode
        :param scores: Mapping from ``reaction_smiles`` to ``(without_fg,
            with_fg)``; must contain ``node.reaction_smiles``.
        :type scores: dict[str, tuple[float, float]]

        :return: Scoring result for this step.
        :rtype: StepResult
        """
        score_without_fg, score_with_fg = scores[node.reaction_smiles]
        return StepResult(
            node=node,
            score_without_fg=score_without_fg,
            score_with_fg=score_with_fg,
        )

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
