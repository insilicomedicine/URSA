from unittest.mock import MagicMock

import pytest

from tests.conftest import _Score
from tests.conftest import make_const_scorer
from tests.conftest import make_scorer
from ursa.basic.result import VariantResult
from ursa.collapsing import PathCollapser
from ursa.configs import PathScoringConfig
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


class TestStepPasses:
    def test_all_pass(self, path_1step, scorer_all_pass):
        result = PathScorer(scorer_all_pass).score(path_1step)
        assert result.all_steps_pass_solv_1 is True
        assert result.all_steps_pass_solv_2 is True
        assert all(sr.passes_solv_1 and sr.passes_solv_2 for sr in result.step_results)

    def test_none_pass(self, path_1step, scorer_all_fail):
        result = PathScorer(scorer_all_fail).score(path_1step)
        assert result.all_steps_pass_solv_1 is False
        assert result.all_steps_pass_solv_2 is False

    def test_error_score_normalised_to_zero(self, path_1step, scorer_error):
        result = PathScorer(scorer_error).score(path_1step)
        assert all(
            sr.score_without_fg == 0.0 and sr.score_with_fg == 0.0
            for sr in result.step_results
        )

    def test_error_score_not_passed(self, path_1step, scorer_error):
        result = PathScorer(scorer_error).score(path_1step)
        assert result.all_steps_pass_solv_1 is False
        assert result.all_steps_pass_solv_2 is False

    def test_fg_split_solv1_passes_solv2_fails(self, path_1step):
        # without_fg > 0 (Solv-1 pass), with_fg == 0 (Solv-2 fail)
        scorer = make_const_scorer((2.0, 0.0))
        result = PathScorer(scorer).score(path_1step)
        assert result.all_steps_pass_solv_1 is True
        assert result.all_steps_pass_solv_2 is False


class TestAggregateMetrics:
    def test_mean_without_fg(self, path_3step_linear):
        scorer = make_scorer(3.0, 1.0, 2.0)
        result = PathScorer(scorer).score(path_3step_linear)
        assert result.mean_score_without_fg == pytest.approx(2.0)
        assert result.mean_score_with_fg == pytest.approx(2.0)

    def test_mean_with_errors_normalised(self, path_3step_linear):
        # -1 -> 0.0, so mean is (3.0 + 0.0 + 2.0) / 3
        scorer = make_scorer(3.0, -1.0, 2.0)
        result = PathScorer(scorer).score(path_3step_linear)
        assert result.mean_score_without_fg == pytest.approx(5.0 / 3)

    def test_distinct_fg_means(self, path_3step_linear):
        scorer = make_scorer((3.0, 1.0), (2.0, 2.0), (1.0, 0.0))
        result = PathScorer(scorer).score(path_3step_linear)
        assert result.mean_score_without_fg == pytest.approx(2.0)
        assert result.mean_score_with_fg == pytest.approx(1.0)


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
        assert len(unique_reactions) < total_steps

        scorer = make_const_scorer(3.0)
        PathScorer(scorer).score_variants(variants)
        assert scorer.evaluate.call_count == len(unique_reactions)

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
            assert b.all_steps_pass_solv_1 == i.all_steps_pass_solv_1
            assert b.all_steps_pass_solv_2 == i.all_steps_pass_solv_2
            assert b.mean_score_without_fg == pytest.approx(i.mean_score_without_fg)
            assert tuple(sr.score_without_fg for sr in b.step_results) == tuple(
                sr.score_without_fg for sr in i.step_results
            )


class TestBuildResultsFromRaw:
    def test_matches_score_variants(self, path_3step_linear):
        variants = PathCollapser().collapse(path_3step_linear)
        unique = _first_seen_reactions(variants)
        raw = {
            reaction: _Score(
                without_functional_groups=float(i + 1),
                with_functional_groups=float(i + 1) / 2,
            )
            for i, reaction in enumerate(unique)
        }

        scorer = MagicMock()
        scorer.evaluate.side_effect = lambda rxn: raw[rxn]
        sequential = PathScorer(scorer).score_variants(variants)

        bulk = PathScorer(MagicMock()).build_results_from_raw(variants, raw)

        for seq, blk in zip(sequential, bulk):
            assert seq.all_steps_pass_solv_1 == blk.all_steps_pass_solv_1
            assert seq.all_steps_pass_solv_2 == blk.all_steps_pass_solv_2
            assert seq.mean_score_with_fg == pytest.approx(blk.mean_score_with_fg)
            assert tuple(s.score_with_fg for s in seq.step_results) == tuple(
                s.score_with_fg for s in blk.step_results
            )

    def test_negative_scores_normalised(self, path_1step):
        variants = (path_1step,)
        reaction = path_1step.get_all_steps()[0].reaction_smiles
        raw = {
            reaction: _Score(
                without_functional_groups=-1.0, with_functional_groups=-1.0
            )
        }
        result = PathScorer(MagicMock()).build_results_from_raw(variants, raw)[0]
        assert result.step_results[0].score_with_fg == 0.0
        assert result.all_steps_pass_solv_2 is False


class TestPassThreshold:
    """A step passes a Solv level iff its score strictly exceeds
    ``PathScoringConfig.pass_threshold`` (configurable, default 0.0)."""

    def test_default_is_zero(self):
        assert PathScoringConfig.pass_threshold.value == 0.0

    def test_score_at_threshold_fails(self, path_1step, monkeypatch):
        # Raise threshold to 1.0: a step scoring exactly 1.0 must NOT pass
        # (strictly-greater rule), even though it passes at the default 0.0.
        assert (
            PathScorer(make_const_scorer(1.0)).score(path_1step).all_steps_pass_solv_2
            is True
        )

        monkeypatch.setattr(PathScoringConfig.pass_threshold, "_value_", 1.0)
        result = PathScorer(make_const_scorer(1.0)).score(path_1step)
        assert result.all_steps_pass_solv_1 is False
        assert result.all_steps_pass_solv_2 is False
        assert result.step_results[0].passes_solv_1 is False

    def test_score_above_threshold_passes(self, path_1step, monkeypatch):
        monkeypatch.setattr(PathScoringConfig.pass_threshold, "_value_", 1.0)
        result = PathScorer(make_const_scorer(2.0)).score(path_1step)
        assert result.all_steps_pass_solv_1 is True
        assert result.all_steps_pass_solv_2 is True
        assert result.step_results[0].passes_solv_2 is True
