from pathlib import Path

import pytest
from chemcensor import ChemCensor

from ursa.basic.node import RetrosyntheticNode
from ursa.basic.path import RetrosyntheticPath
from ursa.basic.result import VariantResult
from ursa.scoring.path_scorer import PathScorer

_FIXTURES_DIR = Path(__file__).parent / "fixtures"
_RC_DB = _FIXTURES_DIR / "rc_db.db"

# Unmapped, canonicalised reaction SMILES extracted from
# chemcensor/tests/unit/extraction/fixtures/reaction_centers_fixtures.json

# fixture 1 — present in rc_db with RC1–RC4 → exact_match (5.0)
_RXN_IN_DB = (
    "COC(=O)c1ccccc1-c1ccc(CN(C(C)=O)[C@@H]2CCCC2N(CC(F)(F)F)C(C)=O)[nH]1"
    ">>"
    "O=C(O)c1ccccc1-c1ccc(CN[C@@H]2CCCC2NCC(F)(F)F)[nH]1"
)

# fixture 6 — not in rc_db → default_reaction_scoring (0.0)
_RXN_NOT_IN_DB = (
    "CC(C)(C)OC(=O)N[C@@H]1CCCC1OS(C)(=O)=O.[N-]=[N+]=[N-]"
    ">>"
    "CC(C)(C)OC(=O)N[C@@H]1CCCC1N=[N+]=[N-]"
)


@pytest.fixture(scope="module")
def chemcensor_scorer():
    if not _RC_DB.exists():
        pytest.skip(f"rc_db.db not found at {_RC_DB}")
    return ChemCensor(db_path=_RC_DB)


def _make_path(rxn_smiles: str, path_id: str = "test") -> RetrosyntheticPath:
    """Build a 1-step path from a forward reaction SMILES ``reactants>>product``."""
    reactants_str, product_str = rxn_smiles.split(">>")
    leaves = tuple(
        RetrosyntheticNode(smiles=s.strip())
        for s in reactants_str.split(".")
        if s.strip()
    )
    root = RetrosyntheticNode(smiles=product_str.strip(), children=leaves)
    return RetrosyntheticPath(path_id=path_id, root=root)


# ── return type ───────────────────────────────────────────────────────────────


class TestReturnType:
    def test_returns_variant_result(self, chemcensor_scorer):
        path = _make_path(_RXN_IN_DB)
        result = PathScorer(chemcensor_scorer).score(path)
        assert isinstance(result, VariantResult)

    def test_path_reference_preserved(self, chemcensor_scorer):
        path = _make_path(_RXN_IN_DB)
        result = PathScorer(chemcensor_scorer).score(path)
        assert result.path is path

    def test_step_count_matches_path(self, chemcensor_scorer):
        path = _make_path(_RXN_IN_DB)
        result = PathScorer(chemcensor_scorer).score(path)
        assert len(result.step_results) == path.num_steps


# ── scoring: reaction in DB ───────────────────────────────────────────────────


class TestScoringInDb:
    def test_exact_match_score(self, chemcensor_scorer):
        path = _make_path(_RXN_IN_DB, path_id="in_db")
        result = PathScorer(chemcensor_scorer).score(path)
        assert result.step_results[0].score_with_fg == pytest.approx(5.0)

    def test_step_passes_solv_2(self, chemcensor_scorer):
        path = _make_path(_RXN_IN_DB, path_id="in_db")
        result = PathScorer(chemcensor_scorer).score(path)
        assert result.step_results[0].passes_solv_2 is True

    def test_all_steps_pass_solv_2(self, chemcensor_scorer):
        path = _make_path(_RXN_IN_DB, path_id="in_db")
        result = PathScorer(chemcensor_scorer).score(path)
        assert result.all_steps_pass_solv_2 is True

    def test_step_passes_solv_1(self, chemcensor_scorer):
        path = _make_path(_RXN_IN_DB, path_id="in_db")
        result = PathScorer(chemcensor_scorer).score(path)
        assert result.step_results[0].passes_solv_1 is True

    def test_mean_with_fg_equals_step_score(self, chemcensor_scorer):
        path = _make_path(_RXN_IN_DB, path_id="in_db")
        result = PathScorer(chemcensor_scorer).score(path)
        assert result.mean_score_with_fg == pytest.approx(
            result.step_results[0].score_with_fg
        )


# ── scoring: reaction not in DB ───────────────────────────────────────────────


class TestScoringNotInDb:
    def test_default_score_zero(self, chemcensor_scorer):
        path = _make_path(_RXN_NOT_IN_DB, path_id="not_in_db")
        result = PathScorer(chemcensor_scorer).score(path)
        assert result.step_results[0].score_with_fg == pytest.approx(0.0)

    def test_step_not_passes_solv_2(self, chemcensor_scorer):
        path = _make_path(_RXN_NOT_IN_DB, path_id="not_in_db")
        result = PathScorer(chemcensor_scorer).score(path)
        assert result.step_results[0].passes_solv_2 is False

    def test_all_steps_pass_solv_2_false(self, chemcensor_scorer):
        path = _make_path(_RXN_NOT_IN_DB, path_id="not_in_db")
        result = PathScorer(chemcensor_scorer).score(path)
        assert result.all_steps_pass_solv_2 is False

    def test_mean_with_fg_is_zero(self, chemcensor_scorer):
        path = _make_path(_RXN_NOT_IN_DB, path_id="not_in_db")
        result = PathScorer(chemcensor_scorer).score(path)
        assert result.mean_score_with_fg == pytest.approx(0.0)


# ── negative score normalisation ─────────────────────────────────────────────


class TestNegativeScoreNormalisation:
    def test_invalid_smiles_normalised_to_zero(self, chemcensor_scorer):
        leaf = RetrosyntheticNode(smiles="C")
        root = RetrosyntheticNode(smiles="CC", children=(leaf,))
        path = RetrosyntheticPath(path_id="invalid", root=root)
        result = PathScorer(chemcensor_scorer).score(path)
        # ChemCensor returns -1.0 for unmappable/invalid reactions
        assert result.step_results[0].score_with_fg >= 0.0
        assert result.step_results[0].score_without_fg >= 0.0

    def test_invalid_smiles_not_passed(self, chemcensor_scorer):
        leaf = RetrosyntheticNode(smiles="C")
        root = RetrosyntheticNode(smiles="CC", children=(leaf,))
        path = RetrosyntheticPath(path_id="invalid", root=root)
        result = PathScorer(chemcensor_scorer).score(path)
        assert result.step_results[0].passes_solv_2 is False


# ── multi-step path ───────────────────────────────────────────────────────────


class TestMultiStep:
    def test_two_steps_both_scored(self, chemcensor_scorer):
        # Build a 2-step path: root -> inter -> [leafA, leafB]
        # Step 1 (inter): reaction not in DB → score 0.0
        # Step 2 (root): reaction in DB → score 5.0
        reactants_str, product_str = _RXN_IN_DB.split(">>")
        leaves = tuple(
            RetrosyntheticNode(smiles=s.strip())
            for s in reactants_str.split(".")
            if s.strip()
        )
        inter = RetrosyntheticNode(smiles=product_str.strip(), children=leaves)
        # Add a wrapping step not in DB
        outer_leaf = RetrosyntheticNode(smiles="C")
        root = RetrosyntheticNode(
            smiles="CC" + product_str.strip(),
            children=(inter, outer_leaf),
        )
        path = RetrosyntheticPath(path_id="two_step", root=root)
        result = PathScorer(chemcensor_scorer).score(path)
        assert len(result.step_results) == 2

    def test_mixed_all_steps_pass_solv_2_false(self, chemcensor_scorer):
        # One step in DB (pass), one not (fail) → all_steps_pass_solv_2 is False
        reactants_str, product_str = _RXN_IN_DB.split(">>")
        leaves = tuple(
            RetrosyntheticNode(smiles=s.strip())
            for s in reactants_str.split(".")
            if s.strip()
        )
        inter = RetrosyntheticNode(smiles=product_str.strip(), children=leaves)
        outer_leaf = RetrosyntheticNode(smiles="C")
        root = RetrosyntheticNode(
            smiles="CC" + product_str.strip(),
            children=(inter, outer_leaf),
        )
        path = RetrosyntheticPath(path_id="mixed", root=root)
        result = PathScorer(chemcensor_scorer).score(path)
        assert result.all_steps_pass_solv_2 is False
