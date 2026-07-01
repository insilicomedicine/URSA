import pytest

from tests.conftest import make_const_scorer
from tests.conftest import make_scorer
from ursa.basic.result import VariantResult
from ursa.collapsing import PathCollapser
from ursa.scoring.path_scorer import PathScorer


def _first_seen_reactions(variants):
    """Reaction smiles in first-seen DFS order (mirrors _score_unique_reactions)."""
    unique: dict[str, None] = {}
    for variant in variants:
        for node in variant.get_all_steps():
            unique.setdefault(node.reaction_smiles, None)
    return unique


@pytest.fixture
def scorer_all_pass():
    return make_const_scorer(3.0)


@pytest.fixture
def scorer_all_fail():
    return make_const_scorer(0.0)


@pytest.fixture
def scorer_error():
    return make_const_scorer(-1.0)


class TestScoreResult:
    def test_returns_variant_result(self, path_1step, scorer_all_pass):
        result = PathScorer(scorer_all_pass).score(path_1step)
        assert isinstance(result, VariantResult)

    def test_step_results_count_matches_steps(self, path_3step_linear, scorer_all_pass):
        result = PathScorer(scorer_all_pass).score(path_3step_linear)
        assert len(result.step_results) == path_3step_linear.num_steps

    def test_path_reference_preserved(self, path_1step, scorer_all_pass):
        result = PathScorer(scorer_all_pass).score(path_1step)
        assert result.path is path_1step


class TestPassedSteps:
    def test_all_passed(self, path_1step, scorer_all_pass):
        result = PathScorer(scorer_all_pass).score(path_1step)
        assert result.all_steps_passed is True
        assert all(sr.passed for sr in result.step_results)

    def test_none_passed(self, path_1step, scorer_all_fail):
        result = PathScorer(scorer_all_fail).score(path_1step)
        assert result.all_steps_passed is False
        assert not any(sr.passed for sr in result.step_results)

    def test_error_score_normalised_to_zero(self, path_1step, scorer_error):
        result = PathScorer(scorer_error).score(path_1step)
        assert all(sr.score == 0.0 for sr in result.step_results)

    def test_error_score_not_passed(self, path_1step, scorer_error):
        result = PathScorer(scorer_error).score(path_1step)
        assert result.all_steps_passed is False


class TestAggregateMetrics:
    def test_percent_found_all_pass(self, path_1step, scorer_all_pass):
        result = PathScorer(scorer_all_pass).score(path_1step)
        assert result.percent_found == pytest.approx(1.0)

    def test_percent_found_none_pass(self, path_1step, scorer_all_fail):
        result = PathScorer(scorer_all_fail).score(path_1step)
        assert result.percent_found == pytest.approx(0.0)

    def test_percent_found_partial(self, path_3step_linear):
        # 2 out of 3 steps pass
        scorer = make_scorer(3.0, 0.0, 2.0)
        result = PathScorer(scorer).score(path_3step_linear)
        assert result.percent_found == pytest.approx(2 / 3)

    def test_chemcensor_per_route_mean(self, path_3step_linear):
        scorer = make_scorer(3.0, 1.0, 2.0)
        result = PathScorer(scorer).score(path_3step_linear)
        assert result.chemcensor_per_route == pytest.approx(2.0)

    def test_chemcensor_per_route_with_errors(self, path_3step_linear):
        # -1 -> 0.0, so mean is (3.0 + 0.0 + 2.0) / 3
        scorer = make_scorer(3.0, -1.0, 2.0)
        result = PathScorer(scorer).score(path_3step_linear)
        assert result.chemcensor_per_route == pytest.approx(5.0 / 3)


class TestScoreVariants:
    def test_one_result_per_variant(self, path_3step_linear, scorer_all_pass):
        variants = PathCollapser().collapse(path_3step_linear)
        results = PathScorer(scorer_all_pass).score_variants(variants)
        assert len(results) == len(variants)
        assert all(isinstance(r, VariantResult) for r in results)
        assert tuple(r.path for r in results) == tuple(variants)

    def test_each_unique_reaction_scored_once(self, path_3step_linear):
        variants = PathCollapser().collapse(path_3step_linear)
        unique_reactions = {
            node.reaction_smiles for v in variants for node in v.get_all_steps()
        }
        total_steps = sum(v.num_steps for v in variants)
        # The collapsed variants share reactions: fewer unique reactions
        # than total steps, otherwise this test proves nothing.
        assert len(unique_reactions) < total_steps

        scorer = make_const_scorer(3.0)
        PathScorer(scorer).score_variants(variants)
        assert scorer.score.call_count == len(unique_reactions)

    def test_matches_independent_scoring(self, path_3step_linear):
        variants = PathCollapser().collapse(path_3step_linear)
        unique_reactions = _first_seen_reactions(variants)
        score_by_reaction = dict(zip(unique_reactions, (3.0, 0.0, 2.0, 0.0, 1.0)))

        batched = PathScorer(
            make_scorer(*(score_by_reaction[r] for r in unique_reactions))
        ).score_variants(variants)

        independent = tuple(
            PathScorer(
                make_scorer(
                    *(score_by_reaction[r] for r in _first_seen_reactions((variant,)))
                )
            ).score(variant)
            for variant in variants
        )

        for b, i in zip(batched, independent):
            assert b.percent_found == pytest.approx(i.percent_found)
            assert b.chemcensor_per_route == pytest.approx(i.chemcensor_per_route)
            assert b.all_steps_passed == i.all_steps_passed
            assert tuple(sr.score for sr in b.step_results) == tuple(
                sr.score for sr in i.step_results
            )
