from typing import Callable

from ..basic import StepResult
from ..basic import VariantResult
from ..configs import PathScoringConfig
from .errors import EmptyVariantsError

ScoreGetter = Callable[[StepResult], float]


def _score_without_fg(step: StepResult) -> float:
    """Solv-1 dimension: functional-group-agnostic step score."""
    return step.score_without_fg


def _score_with_fg(step: StepResult) -> float:
    """Solv-2 dimension: functional-group-aware step score."""
    return step.score_with_fg


class BestPathSelector:
    """Selects the best collapsed variant for a given Solv dimension.

    Selection is parametrised by a ``score_getter`` that extracts the
    relevant per-step score (functional-group-agnostic for Solv-1,
    functional-group-aware for Solv-2). The three-level priority is:

    1. **Minimum failed steps** — fewest steps scoring at or below
       :attr:`PathScoringConfig.pass_threshold` in the chosen dimension.
    2. **Maximum average score** — highest mean score in that dimension.
    3. **Shortest path** — fewest steps, as the most compressed synthesis.

    Because selection minimises failed steps first, a route satisfies the
    corresponding Solv level *iff* the selected variant has zero failed
    steps — i.e. selection doubles as the existence check over variants.
    """

    def select(
        self,
        variants: tuple[VariantResult, ...],
        score_getter: ScoreGetter,
    ) -> VariantResult:
        """Return the best variant under ``score_getter``.

        :param variants: All scored variants of a retrosynthetic path
            (original and collapsed forms). Must not be empty.
        :type variants: tuple[VariantResult, ...]
        :param score_getter: Extracts the per-step score to rank by.
        :type score_getter: ScoreGetter

        :return: The best variant for the chosen dimension.
        :rtype: VariantResult
        """
        if not variants:
            raise EmptyVariantsError("unknown")
        return min(
            variants,
            key=lambda v: (
                self._count_failed_steps(v, score_getter),
                -self._average_score(v, score_getter),
                v.path.num_steps,
            ),
        )

    def select_for_solv_1(self, variants: tuple[VariantResult, ...]) -> VariantResult:
        """Select the best variant for Solv-1 (functional-group-agnostic)."""
        return self.select(variants, _score_without_fg)

    def select_for_solv_2(self, variants: tuple[VariantResult, ...]) -> VariantResult:
        """Select the best variant for Solv-2 (functional-group-aware)."""
        return self.select(variants, _score_with_fg)

    def _count_failed_steps(
        self, variant: VariantResult, score_getter: ScoreGetter
    ) -> int:
        """Return the number of failed steps in ``variant``.

        A step fails when its score does not strictly exceed
        :attr:`PathScoringConfig.pass_threshold`.

        :param variant: The variant to evaluate.
        :type variant: VariantResult
        :param score_getter: Extracts the per-step score to test.
        :type score_getter: ScoreGetter

        :return: Count of failed steps (score at or below ``pass_threshold``)
            in the chosen dimension.
        :rtype: int
        """
        return sum(
            1
            for sr in variant.step_results
            if score_getter(sr) <= PathScoringConfig.pass_threshold.value
        )

    def _average_score(
        self, variant: VariantResult, score_getter: ScoreGetter
    ) -> float:
        """Return the mean step score of ``variant`` in the chosen dimension.

        :param variant: The variant to evaluate.
        :type variant: VariantResult
        :param score_getter: Extracts the per-step score to average.
        :type score_getter: ScoreGetter

        :return: Mean score; ``0.0`` for an empty variant.
        :rtype: float
        """
        if not variant.step_results:
            return 0.0
        return sum(score_getter(sr) for sr in variant.step_results) / len(
            variant.step_results
        )
