import pytest

from tests.conftest import make_variant
from ursa.basic.building_block import BuildingBlock
from ursa.basic.result import DatasetMetrics
from ursa.basic.result import PathResult
from ursa.metrics.calculator import DatasetMetricsCalculator


def make_path_result(
    path,
    bb_found: bool,
    scores: tuple[float, ...],
    is_consistent: bool = True,
) -> PathResult:
    bb = (BuildingBlock(smiles="CC", found_in_catalog=bb_found),)
    all_bb_found = bb_found
    best_variant = make_variant(path, scores)
    is_route_solved = is_consistent and all_bb_found and best_variant.all_steps_passed
    return PathResult(
        original_path=path,
        is_consistent=is_consistent,
        starting_materials=bb,
        all_bb_found=all_bb_found,
        best_variant=best_variant,
        is_route_solved=is_route_solved,
    )


@pytest.fixture
def calculator() -> DatasetMetricsCalculator:
    return DatasetMetricsCalculator()


class TestCalculate:
    def test_returns_dataset_metrics(self, calculator, path_1step):
        pr = make_path_result(path_1step, True, (3.0,))
        result = calculator.calculate((pr,), total_molecules=1)
        assert isinstance(result, DatasetMetrics)

    def test_total_molecules_from_arg(self, calculator, path_1step):
        pr = make_path_result(path_1step, True, (3.0,))
        result = calculator.calculate((pr,), total_molecules=10)
        assert result.total_molecules == 10

    def test_molecules_with_route(self, calculator, path_1step):
        pr = make_path_result(path_1step, True, (3.0,))
        result = calculator.calculate((pr, pr), total_molecules=5)
        assert result.molecules_with_route == 2

    def test_solved_routes_count(self, calculator, path_1step):
        solved = make_path_result(path_1step, True, (3.0,))
        unsolved = make_path_result(path_1step, True, (0.0,))
        result = calculator.calculate((solved, unsolved), total_molecules=2)
        assert result.solved_routes == 1

    def test_solv_2(self, calculator, path_1step):
        solved = make_path_result(path_1step, True, (3.0,))
        result = calculator.calculate((solved,), total_molecules=4)
        assert result.solv_2 == pytest.approx(1 / 4)

    def test_solv_2_zero_molecules(self, calculator):
        result = calculator.calculate((), total_molecules=0)
        assert result.solv_2 == pytest.approx(0.0)

    def test_passed_steps(self, calculator, path_3step_linear):
        pr = make_path_result(path_3step_linear, True, (3.0, 0.0, 2.0))
        result = calculator.calculate((pr,), total_molecules=1)
        assert result.passed_steps == 2

    def test_total_steps(self, calculator, path_3step_linear):
        pr = make_path_result(path_3step_linear, True, (3.0, 1.0, 2.0))
        result = calculator.calculate((pr,), total_molecules=1)
        assert result.total_steps == 3

    def test_mean_chemcensor_score(self, calculator, path_3step_linear):
        pr = make_path_result(path_3step_linear, True, (3.0, 1.0, 2.0))
        result = calculator.calculate((pr,), total_molecules=1)
        assert result.mean_chemcensor_score == pytest.approx(2.0)

    def test_mean_chemcensor_score_no_steps(self, calculator):
        result = calculator.calculate((), total_molecules=1)
        assert result.mean_chemcensor_score == pytest.approx(0.0)

    def test_multiple_paths_aggregate(self, calculator, path_1step, path_3step_linear):
        pr1 = make_path_result(path_1step, True, (3.0,))
        pr2 = make_path_result(path_3step_linear, True, (2.0, 0.0, 1.0))
        result = calculator.calculate((pr1, pr2), total_molecules=2)
        assert result.total_steps == 4
        assert result.passed_steps == 3  # 1 + 2
        assert result.solved_routes == 1  # only pr1 (pr2 has a failed step)
