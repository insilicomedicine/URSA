from ..basic import VariantResult
from ..configs import PathScoringConfig
from .errors import EmptyVariantsError


class BestPathSelector:
    """Selects the best-scoring variant from a set of collapsed variants.

    The selection follows a three-level priority:

    1. **Minimum failed steps** — variants with the fewest steps where
       ``score == 0`` (not found in the database) are preferred.
    2. **Maximum average score** — among equally failing variants, the
       one with the highest ``chemcensor_per_route`` is preferred.
    3. **Shortest path** — among variants still tied, the shortest
       (fewest steps) is preferred, as it represents the most compressed
       description of the synthesis.
    """

    def select(self, variants: tuple[VariantResult, ...]) -> VariantResult:
        """Return the best variant from ``variants``.

        :param variants: All scored variants of a retrosynthetic path,
            including the original and all collapsed forms. Must not be
            empty.
        :type variants: tuple[VariantResult, ...]

        :return: The variant with the best score according to the
            three-level priority rule.
        :rtype: VariantResult
        """
        if not variants:
            raise EmptyVariantsError("unknown")
        return min(
            variants,
            key=lambda v: (
                self._count_failed_steps(v),
                -self._average_score(v),
                v.path.num_steps,
            ),
        )

    def _count_failed_steps(self, variant: VariantResult) -> int:
        """Return the number of steps in ``variant`` with ``score == 0``.

        Steps scored as SIS or tautomerization are excluded (they are
        treated as automatically passing with score 1.0).

        :param variant: The variant to evaluate.
        :type variant: VariantResult

        :return: Count of failed (score == 0) steps.
        :rtype: int
        """
        return sum(
            1
            for sr in variant.step_results
            if sr.score == PathScoringConfig.failed_step_score.value
        )

    def _average_score(self, variant: VariantResult) -> float:
        """Return the mean reaction score across all steps of ``variant``.

        Empty variants (no steps) return 0.0.

        :param variant: The variant to evaluate.
        :type variant: VariantResult

        :return: Mean score in ``[0.0, max_scorer_value]``.
        :rtype: float
        """
        if not variant.step_results:
            return 0.0
        return sum(sr.score for sr in variant.step_results) / len(variant.step_results)
