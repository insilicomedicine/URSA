import pytest

from tests.conftest import make_variant
from ursa.scoring.best_path_selector import BestPathSelector
from ursa.scoring.errors import EmptyVariantsError


@pytest.fixture
def selector() -> BestPathSelector:
    return BestPathSelector()


class TestSelect:
    def test_single_variant_returned(self, selector, path_1step):
        vr = make_variant(path_1step, (3.0,))
        assert selector.select((vr,)) is vr

    def test_empty_variants_raises(self, selector):
        with pytest.raises(EmptyVariantsError):
            selector.select(())

    def test_prefers_fewer_failed_steps(self, selector, path_3step_linear):
        # vr0: 0 failed (scores 3,2,1); vr1: 1 failed (scores 3,0,2)
        vr0 = make_variant(path_3step_linear, (3.0, 2.0, 1.0))
        vr1 = make_variant(path_3step_linear, (3.0, 0.0, 2.0))
        assert selector.select((vr1, vr0)) is vr0

    def test_tiebreak_by_higher_average(self, selector, path_3step_linear):
        # Both have 0 failed; vr_high has higher avg
        vr_low = make_variant(path_3step_linear, (1.0, 1.0, 1.0))
        vr_high = make_variant(path_3step_linear, (3.0, 2.0, 1.0))
        assert selector.select((vr_low, vr_high)) is vr_high

    def test_tiebreak_by_shorter_path(self, selector, path_3step_linear, path_1step):
        # Both 0 failed, same avg; shorter path wins
        vr_long = make_variant(path_3step_linear, (2.0, 2.0, 2.0))
        vr_short = make_variant(path_1step, (2.0,))
        assert selector.select((vr_long, vr_short)) is vr_short


class TestCountFailedSteps:
    def test_no_failed(self, selector, path_1step):
        vr = make_variant(path_1step, (3.0,))
        assert selector._count_failed_steps(vr) == 0

    def test_all_failed(self, selector, path_3step_linear):
        vr = make_variant(path_3step_linear, (0.0, 0.0, 0.0))
        assert selector._count_failed_steps(vr) == 3

    def test_partial_failed(self, selector, path_3step_linear):
        vr = make_variant(path_3step_linear, (3.0, 0.0, 2.0))
        assert selector._count_failed_steps(vr) == 1

    def test_sis_score_not_counted_as_failed(self, selector, path_1step):
        # SIS stored as 1.0 (>= pass_threshold, score != 0)
        vr = make_variant(path_1step, (1.0,))
        assert selector._count_failed_steps(vr) == 0


class TestAverageScore:
    def test_empty_variant(self, selector, path_1step):
        from ursa.basic.result import VariantResult

        vr = VariantResult(
            path=path_1step,
            step_results=(),
            percent_found=0.0,
            chemcensor_per_route=0.0,
            all_steps_passed=True,
        )
        assert selector._average_score(vr) == pytest.approx(0.0)

    def test_single_step(self, selector, path_1step):
        vr = make_variant(path_1step, (4.0,))
        assert selector._average_score(vr) == pytest.approx(4.0)

    def test_multiple_steps_mean(self, selector, path_3step_linear):
        vr = make_variant(path_3step_linear, (3.0, 1.0, 2.0))
        assert selector._average_score(vr) == pytest.approx(2.0)
