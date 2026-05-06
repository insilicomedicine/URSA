from ..basic import RetrosyntheticNode
from ..basic import RetrosyntheticPath
from ..basic import StepResult
from ..basic import VariantResult
from ..configs import PathScoringConfig
from .base import ReactionScorer


class PathScorer:
    """Scores all retrosynthetic steps of a single path variant.

    Traverses the tree in DFS pre-order, calls the injected
    :class:`~ursa.ReactionScorer` on each internal node's
    ``reaction_smiles``, and aggregates the per-step results into a
    :class:`~ursa.VariantResult`.

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

        Traverses the tree in DFS pre-order. For each internal node,
        calls :meth:`_score_node`. Aggregates ``percent_found``,
        ``chemcensor_per_route``, and
        ``all_steps_passed`` from the individual step results.

        :param path: The retrosynthetic path variant to score.
        :type path: RetrosyntheticPath

        :return: Aggregated scoring result for the variant.
        :rtype: VariantResult
        """
        step_results = tuple(self._score_node(n) for n in path.get_all_steps())
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

    def _score_node(self, node: RetrosyntheticNode) -> StepResult:
        """Score a single internal node.

        Calls the injected :class:`~ursa.ReactionScorer` with
        ``node.reaction_smiles`` and constructs a :class:`~ursa.StepResult`.
        A step is considered *passed* when its score is at or above
        :attr:`PathScoringConfig.pass_threshold`.

        :param node: Internal node whose ``reaction_smiles`` is scored.
        :type node: RetrosyntheticNode

        :return: Scoring result for this step.
        :rtype: StepResult
        """
        raw = self._scorer.score(node.reaction_smiles)
        score = PathScoringConfig.failed_step_score.value if raw < 0 else float(raw)
        passed = score >= PathScoringConfig.pass_threshold.value
        return StepResult(node=node, score=score, passed=passed)
