import pytest

from tests.conftest import make_path_result
from ursa.basic.result import DatasetMetrics
from ursa.metrics.calculator import DatasetMetricsCalculator


@pytest.fixture
def calculator() -> DatasetMetricsCalculator:
    return DatasetMetricsCalculator()


class TestCalculate:
    def test_returns_dataset_metrics(self, calculator, path_1step):
        pr = make_path_result(path_1step, (3.0,))
        result = calculator.calculate((pr,), total_molecules=1)
        assert isinstance(result, DatasetMetrics)

    def test_total_molecules_from_arg(self, calculator, path_1step):
        pr = make_path_result(path_1step, (3.0,))
        result = calculator.calculate((pr,), total_molecules=10)
        assert result.total_molecules == 10

    def test_molecules_with_route(self, calculator, path_1step):
        pr = make_path_result(path_1step, (3.0,))
        result = calculator.calculate((pr, pr), total_molecules=5)
        assert result.molecules_with_route == 2

    def test_routes_solv_2_count(self, calculator, path_1step):
        passing = make_path_result(path_1step, (3.0,))
        failing = make_path_result(path_1step, (0.0,))
        result = calculator.calculate((passing, failing), total_molecules=2)
        assert result.routes_solv_2 == 1

    def test_solv_2_rate(self, calculator, path_1step):
        passing = make_path_result(path_1step, (3.0,))
        result = calculator.calculate((passing,), total_molecules=4)
        assert result.solv_2 == pytest.approx(1 / 4)

    def test_solv_rates_zero_molecules(self, calculator):
        result = calculator.calculate((), total_molecules=0)
        assert result.solv_0 == pytest.approx(0.0)
        assert result.solv_1 == pytest.approx(0.0)
        assert result.solv_2 == pytest.approx(0.0)

    def test_solv_hierarchy_counts(self, calculator, path_3step_linear):
        # Solv-1 passes (without-fg all > 0) but Solv-2 fails (a with-fg zero).
        pr = make_path_result(
            path_3step_linear,
            scores=(3.0, 2.0, 1.0),
            scores_with_fg=(3.0, 0.0, 1.0),
        )
        result = calculator.calculate((pr,), total_molecules=1)
        assert result.routes_solv_0 == 1
        assert result.routes_solv_1 == 1
        assert result.routes_solv_2 == 0

    def test_solv_0_requires_bb_and_consistency(self, calculator, path_1step):
        bb_missing = make_path_result(path_1step, (3.0,), bb_found=False)
        result = calculator.calculate((bb_missing,), total_molecules=1)
        assert result.routes_solv_0 == 0
        assert result.routes_solv_1 == 0
        assert result.routes_solv_2 == 0

    def test_mean_scores(self, calculator, path_3step_linear):
        pr = make_path_result(
            path_3step_linear,
            scores=(3.0, 1.0, 2.0),
            scores_with_fg=(2.0, 2.0, 2.0),
        )
        result = calculator.calculate((pr,), total_molecules=1)
        assert result.mean_score_without_fg == pytest.approx(2.0)
        assert result.mean_score_with_fg == pytest.approx(2.0)

    def test_mean_scores_no_steps(self, calculator):
        result = calculator.calculate((), total_molecules=1)
        assert result.mean_score_without_fg == pytest.approx(0.0)
        assert result.mean_score_with_fg == pytest.approx(0.0)

    def test_multiple_paths_aggregate(self, calculator, path_1step, path_3step_linear):
        pr1 = make_path_result(path_1step, (3.0,))
        pr2 = make_path_result(path_3step_linear, (2.0, 0.0, 1.0))
        result = calculator.calculate((pr1, pr2), total_molecules=2)
        assert result.molecules_with_route == 2
        assert result.routes_solv_2 == 1  # only pr1 (pr2 has a failed step)
