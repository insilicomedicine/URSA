import pytest

from tests.conftest import make_variant
from ursa.configs import PathScoringConfig
from ursa.scoring.best_path_selector import _score_with_fg
from ursa.scoring.best_path_selector import BestPathSelector
from ursa.scoring.errors import EmptyVariantsError


@pytest.fixture
def selector() -> BestPathSelector:
    return BestPathSelector()


class TestSelect:
    def test_single_variant_returned(self, selector, path_1step):
        vr = make_variant(path_1step, (3.0,))
        assert selector.select_for_solv_2((vr,)) is vr

    def test_empty_variants_raises(self, selector):
        with pytest.raises(EmptyVariantsError):
            selector.select_for_solv_2(())

    def test_prefers_fewer_failed_steps(self, selector, path_3step_linear):
        vr0 = make_variant(path_3step_linear, (3.0, 2.0, 1.0))
        vr1 = make_variant(path_3step_linear, (3.0, 0.0, 2.0))
        assert selector.select_for_solv_2((vr1, vr0)) is vr0

    def test_tiebreak_by_higher_average(self, selector, path_3step_linear):
        vr_low = make_variant(path_3step_linear, (1.0, 1.0, 1.0))
        vr_high = make_variant(path_3step_linear, (3.0, 2.0, 1.0))
        assert selector.select_for_solv_2((vr_low, vr_high)) is vr_high

    def test_tiebreak_by_shorter_path(self, selector, path_3step_linear, path_1step):
        vr_long = make_variant(path_3step_linear, (2.0, 2.0, 2.0))
        vr_short = make_variant(path_1step, (2.0,))
        assert selector.select_for_solv_2((vr_long, vr_short)) is vr_short

    def test_solv_1_and_solv_2_can_pick_different_variants(
        self, selector, path_3step_linear
    ):
        # vr_a passes Solv-1 fully (no without-fg zeros) but fails Solv-2;
        # vr_b is the reverse-ish. Selection should pick each independently.
        vr_a = make_variant(
            path_3step_linear,
            scores=(2.0, 2.0, 2.0),
            scores_with_fg=(2.0, 0.0, 2.0),
        )
        vr_b = make_variant(
            path_3step_linear,
            scores=(2.0, 0.0, 2.0),
            scores_with_fg=(2.0, 2.0, 2.0),
        )
        assert selector.select_for_solv_1((vr_a, vr_b)) is vr_a
        assert selector.select_for_solv_2((vr_a, vr_b)) is vr_b


class TestCountFailedSteps:
    def test_no_failed(self, selector, path_1step):
        vr = make_variant(path_1step, (3.0,))
        assert selector._count_failed_steps(vr, _score_with_fg) == 0

    def test_all_failed(self, selector, path_3step_linear):
        vr = make_variant(path_3step_linear, (0.0, 0.0, 0.0))
        assert selector._count_failed_steps(vr, _score_with_fg) == 3

    def test_partial_failed(self, selector, path_3step_linear):
        vr = make_variant(path_3step_linear, (3.0, 0.0, 2.0))
        assert selector._count_failed_steps(vr, _score_with_fg) == 1

    def test_positive_score_not_counted_as_failed(self, selector, path_1step):
        vr = make_variant(path_1step, (1.0,))
        assert selector._count_failed_steps(vr, _score_with_fg) == 0

    def test_respects_pass_threshold(self, selector, path_3step_linear, monkeypatch):
        # At threshold 1.5, scores <= 1.5 count as failed (here: 1.0 and 1.5),
        # while 2.0 passes.
        vr = make_variant(path_3step_linear, (1.0, 1.5, 2.0))
        assert selector._count_failed_steps(vr, _score_with_fg) == 0

        monkeypatch.setattr(PathScoringConfig.pass_threshold, "_value_", 1.5)
        assert selector._count_failed_steps(vr, _score_with_fg) == 2


class TestAverageScore:
    def test_empty_variant(self, selector, path_1step):
        vr = make_variant(path_1step, ())
        assert selector._average_score(vr, _score_with_fg) == pytest.approx(0.0)

    def test_single_step(self, selector, path_1step):
        vr = make_variant(path_1step, (4.0,))
        assert selector._average_score(vr, _score_with_fg) == pytest.approx(4.0)

    def test_multiple_steps_mean(self, selector, path_3step_linear):
        vr = make_variant(path_3step_linear, (3.0, 1.0, 2.0))
        assert selector._average_score(vr, _score_with_fg) == pytest.approx(2.0)
