import pytest

from tests.conftest import make_const_scorer
from tests.conftest import make_scorer
from ursa.basic.result import VariantResult
from ursa.scoring.path_scorer import PathScorer


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
