from pathlib import Path

import pytest

from chemcensor import ChemCensor
from ursa import DatasetResult
from ursa import PathResult
from ursa import Ursa
from ursa.basic.node import RetrosyntheticNode
from ursa.basic.path import RetrosyntheticPath

_FIXTURES_DIR = Path(__file__).parent / "fixtures"
_RC_DB = _FIXTURES_DIR / "rc_db.db"

# fixture 1 from chemcensor fixtures — present in rc_db with RC1–RC4
_RXN_IN_DB = (
    "COC(=O)c1ccccc1-c1ccc(CN(C(C)=O)[C@@H]2CCCC2N(CC(F)(F)F)C(C)=O)[nH]1"
    ">>"
    "O=C(O)c1ccccc1-c1ccc(CN[C@@H]2CCCC2NCC(F)(F)F)[nH]1"
)

# fixture 6 — not in rc_db → score 0.0
_RXN_NOT_IN_DB = (
    "CC(C)(C)OC(=O)N[C@@H]1CCCC1OS(C)(=O)=O.[N-]=[N+]=[N-]"
    ">>"
    "CC(C)(C)OC(=O)N[C@@H]1CCCC1N=[N+]=[N-]"
)


# ── helpers ───────────────────────────────────────────────────────────────────


def _make_path(rxn_smiles: str, path_id: str = "test") -> RetrosyntheticPath:
    """Build a minimal 1-step path from ``reactants>>product``."""
    reactants_str, product_str = rxn_smiles.split(">>")
    leaves = tuple(
        RetrosyntheticNode(smiles=s.strip())
        for s in reactants_str.split(".")
        if s.strip()
    )
    root = RetrosyntheticNode(smiles=product_str.strip(), children=leaves)
    return RetrosyntheticPath(path_id=path_id, root=root)


def _catalog_file(tmp_path: Path, smiles_lines: list[str]) -> Path:
    p = tmp_path / "catalog.smi"
    p.write_text("\n".join(smiles_lines) + "\n")
    return p


# ── fixtures ──────────────────────────────────────────────────────────────────


@pytest.fixture(scope="module")
def scorer():
    if not _RC_DB.exists():
        pytest.skip(f"rc_db.db not found at {_RC_DB}")
    return ChemCensor(db_path=_RC_DB)


@pytest.fixture
def path_in_db():
    return _make_path(_RXN_IN_DB, path_id="in_db")


@pytest.fixture
def path_not_in_db():
    return _make_path(_RXN_NOT_IN_DB, path_id="not_in_db")


@pytest.fixture
def ursa_all_bb_found(scorer, tmp_path, path_in_db):
    """Ursa with a catalog that contains all leaves of path_in_db."""
    leaves = [n.canonical_smiles for n in path_in_db.get_starting_materials()]
    catalog = _catalog_file(tmp_path, leaves)
    return Ursa(scorer=scorer, bb_catalog_path=catalog)


@pytest.fixture
def ursa_no_bb_found(scorer, tmp_path):
    """Ursa with an empty catalog."""
    catalog = _catalog_file(tmp_path, [])
    return Ursa(scorer=scorer, bb_catalog_path=catalog)


# ── Ursa.score: return type ───────────────────────────────────────────────────


class TestScoreReturnType:
    def test_returns_path_result(self, ursa_all_bb_found, path_in_db):
        result = ursa_all_bb_found.score(path_in_db)
        assert isinstance(result, PathResult)

    def test_original_path_preserved(self, ursa_all_bb_found, path_in_db):
        result = ursa_all_bb_found.score(path_in_db)
        assert result.original_path is path_in_db


# ── Ursa.score: consistency ───────────────────────────────────────────────────


class TestScoreConsistency:
    def test_valid_path_is_consistent(self, ursa_all_bb_found, path_in_db):
        result = ursa_all_bb_found.score(path_in_db)
        assert result.is_consistent is True

    def test_invalid_smiles_path_is_not_consistent(self, scorer, tmp_path):
        catalog = _catalog_file(tmp_path, [])
        ursa = Ursa(scorer=scorer, bb_catalog_path=catalog)
        invalid_leaf = RetrosyntheticNode(smiles="NOT_VALID!!!")
        root = RetrosyntheticNode(
            smiles="CCO", children=(invalid_leaf, RetrosyntheticNode(smiles="O"))
        )
        path = RetrosyntheticPath(path_id="bad", root=root)
        result = ursa.score(path)
        assert result.is_consistent is False


# ── Ursa.score: building blocks ───────────────────────────────────────────────


class TestScoreBuildingBlocks:
    def test_all_bb_found_when_catalog_complete(self, ursa_all_bb_found, path_in_db):
        result = ursa_all_bb_found.score(path_in_db)
        assert result.all_bb_found is True

    def test_no_bb_found_with_empty_catalog(self, ursa_no_bb_found, path_in_db):
        result = ursa_no_bb_found.score(path_in_db)
        assert result.all_bb_found is False

    def test_starting_materials_count(self, ursa_all_bb_found, path_in_db):
        result = ursa_all_bb_found.score(path_in_db)
        expected = len(path_in_db.get_starting_materials())
        assert len(result.starting_materials) == expected


# ── Ursa.score: is_route_solved ───────────────────────────────────────────────


class TestScoreIsRouteSolved:
    def test_solved_when_consistent_bb_found_all_steps_pass(
        self, ursa_all_bb_found, path_in_db
    ):
        result = ursa_all_bb_found.score(path_in_db)
        # path_in_db uses _RXN_IN_DB (exact match, score 5.0 → passed)
        assert result.is_route_solved is True

    def test_not_solved_when_bb_missing(self, ursa_no_bb_found, path_in_db):
        result = ursa_no_bb_found.score(path_in_db)
        assert result.is_route_solved is False

    def test_not_solved_when_step_fails(self, scorer, tmp_path, path_not_in_db):
        leaves = [n.canonical_smiles for n in path_not_in_db.get_starting_materials()]
        catalog = _catalog_file(tmp_path, leaves)
        ursa = Ursa(scorer=scorer, bb_catalog_path=catalog)
        result = ursa.score(path_not_in_db)
        # _RXN_NOT_IN_DB scores 0.0 → step fails → not solved
        assert result.is_route_solved is False


# ── Ursa.score: best_variant ──────────────────────────────────────────────────


class TestScoreBestVariant:
    def test_best_variant_score_in_db(self, ursa_all_bb_found, path_in_db):
        result = ursa_all_bb_found.score(path_in_db)
        assert result.best_variant.step_results[0].score == pytest.approx(5.0)

    def test_best_variant_score_not_in_db(self, ursa_no_bb_found, path_not_in_db):
        result = ursa_no_bb_found.score(path_not_in_db)
        assert result.best_variant.step_results[0].score == pytest.approx(0.0)


# ── Ursa.score_dataset ────────────────────────────────────────────────────────


class TestScoreDataset:
    def test_returns_dataset_result(self, ursa_all_bb_found, path_in_db):
        result = ursa_all_bb_found.score_dataset([path_in_db])
        assert isinstance(result, DatasetResult)

    def test_path_results_count(self, ursa_all_bb_found, path_in_db, path_not_in_db):
        result = ursa_all_bb_found.score_dataset([path_in_db, path_not_in_db])
        assert len(result.path_results) == 2

    def test_total_molecules_defaults_to_len(self, ursa_all_bb_found, path_in_db):
        result = ursa_all_bb_found.score_dataset([path_in_db, path_in_db])
        assert result.metrics.total_molecules == 2

    def test_total_molecules_explicit(self, ursa_all_bb_found, path_in_db):
        result = ursa_all_bb_found.score_dataset([path_in_db], total_molecules=10)
        assert result.metrics.total_molecules == 10

    def test_solved_routes_count(self, scorer, tmp_path, path_in_db, path_not_in_db):
        leaves = [n.canonical_smiles for n in path_in_db.get_starting_materials()]
        catalog = _catalog_file(tmp_path, leaves)
        ursa = Ursa(scorer=scorer, bb_catalog_path=catalog)
        result = ursa.score_dataset([path_in_db, path_not_in_db], total_molecules=2)
        # path_in_db: solved; path_not_in_db: not solved (step fails)
        assert result.metrics.solved_routes == 1

    def test_solv_2(self, scorer, tmp_path, path_in_db, path_not_in_db):
        leaves = [n.canonical_smiles for n in path_in_db.get_starting_materials()]
        catalog = _catalog_file(tmp_path, leaves)
        ursa = Ursa(scorer=scorer, bb_catalog_path=catalog)
        result = ursa.score_dataset([path_in_db, path_not_in_db], total_molecules=4)
        assert result.metrics.solv_2 == pytest.approx(1 / 4)

    def test_empty_dataset(self, ursa_all_bb_found):
        result = ursa_all_bb_found.score_dataset([])
        assert result.metrics.total_molecules == 0
        assert result.metrics.solved_routes == 0
        assert result.metrics.solv_2 == pytest.approx(0.0)


# ── Ursa.score_dataset with target_smiles ─────────────────────────────────────


class TestScoreDatasetWithTargets:
    def test_off_target_path_discarded(self, scorer, tmp_path, path_in_db):
        """A path whose root is not in target_smiles is ignored."""
        other_path = _make_path(_RXN_NOT_IN_DB, path_id="other")
        leaves = [n.canonical_smiles for n in path_in_db.get_starting_materials()]
        catalog = _catalog_file(tmp_path, leaves)
        ursa = Ursa(scorer=scorer, bb_catalog_path=catalog)

        product_in_db = _RXN_IN_DB.split(">>")[1].strip()
        result = ursa.score_dataset(
            [path_in_db, other_path],
            target_smiles=[product_in_db],
        )
        assert len(result.path_results) == 1

    def test_total_molecules_derived_from_targets(self, scorer, tmp_path, path_in_db):
        """total_molecules equals the number of unique canonical target SMILES."""
        leaves = [n.canonical_smiles for n in path_in_db.get_starting_materials()]
        catalog = _catalog_file(tmp_path, leaves)
        ursa = Ursa(scorer=scorer, bb_catalog_path=catalog)

        product_in_db = _RXN_IN_DB.split(">>")[1].strip()
        result = ursa.score_dataset(
            [path_in_db],
            target_smiles=[product_in_db, "CCCO"],
        )
        assert result.metrics.total_molecules == 2

    def test_target_without_path_lowers_solvability(self, scorer, tmp_path, path_in_db):
        """A target with no matching path counts in denominator but not numerator."""
        leaves = [n.canonical_smiles for n in path_in_db.get_starting_materials()]
        catalog = _catalog_file(tmp_path, leaves)
        ursa = Ursa(scorer=scorer, bb_catalog_path=catalog)

        product_in_db = _RXN_IN_DB.split(">>")[1].strip()
        result = ursa.score_dataset(
            [path_in_db],
            target_smiles=[product_in_db, "CCCO"],
        )
        # 1 solved out of 2 targets
        assert result.metrics.solv_2 == pytest.approx(result.metrics.solved_routes / 2)

    def test_multiple_paths_same_target_best_kept(
        self, scorer, tmp_path, path_in_db, path_not_in_db
    ):
        """When two paths share a target root, only the best PathResult is kept."""
        product_in_db = _RXN_IN_DB.split(">>")[1].strip()
        # Rebuild path_not_in_db with the same root as path_in_db
        #   so both paths have the same canonical root SMILES
        leaves_in = [n.canonical_smiles for n in path_in_db.get_starting_materials()]
        catalog = _catalog_file(tmp_path, leaves_in)
        ursa = Ursa(scorer=scorer, bb_catalog_path=catalog)

        # Make a "worse" path with same root, but wrong leaves → step fails
        wrong_leaf = RetrosyntheticNode(smiles="C")
        bad_root = RetrosyntheticNode(
            smiles=product_in_db,
            children=(wrong_leaf,),
        )
        bad_path = RetrosyntheticPath(path_id="bad", root=bad_root)

        result = ursa.score_dataset(
            [path_in_db, bad_path],
            target_smiles=[product_in_db],
        )
        assert len(result.path_results) == 1
        # The good path should have been selected (exact match score)
        assert result.path_results[0].best_variant.step_results[
            0
        ].score == pytest.approx(5.0)
