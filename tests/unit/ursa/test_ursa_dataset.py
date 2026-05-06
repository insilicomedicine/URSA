from pathlib import Path

import pytest
from rdkit import Chem

from tests.conftest import make_const_scorer
from tests.conftest import make_variant
from ursa.basic.building_block import BuildingBlock
from ursa.basic.node import RetrosyntheticNode
from ursa.basic.path import RetrosyntheticPath
from ursa.basic.result import DatasetResult
from ursa.basic.result import PathResult
from ursa.ursa import Ursa


# ── helpers ───────────────────────────────────────────────────────────────────


def _make_ursa(
    tmp_path: Path,
    score: float = 3.0,
    catalog_lines: list[str] | None = None,
) -> Ursa:
    """Build an Ursa instance backed by a const-scorer mock.

    The catalog file is created before Ursa is initialised so that
    BuildingBlockChecker can open it immediately.
    """
    cat = tmp_path / "catalog.smi"
    cat.write_text("\n".join(catalog_lines or []) + "\n")
    return Ursa(scorer=make_const_scorer(score), bb_catalog_path=str(cat))


def _make_path(
    root_smiles: str, leaf_smiles: str = "CC", path_id: str = "p"
) -> RetrosyntheticPath:
    root = RetrosyntheticNode(
        smiles=root_smiles,
        children=(RetrosyntheticNode(smiles=leaf_smiles),),
    )
    return RetrosyntheticPath(path_id=path_id, root=root)


def _make_path_result(
    path: RetrosyntheticPath,
    scores: tuple[float, ...],
    is_route_solved: bool = True,
) -> PathResult:
    vr = make_variant(path, scores)
    return PathResult(
        original_path=path,
        is_consistent=True,
        starting_materials=(BuildingBlock(smiles="CC", found_in_catalog=True),),
        all_bb_found=True,
        best_variant=vr,
        is_route_solved=is_route_solved,
    )


def _canonical(smiles: str) -> str:
    mol = Chem.MolFromSmiles(smiles)
    return Chem.MolToSmiles(mol) if mol is not None else ""


# ── GroupPathsByTarget ────────────────────────────────────────────────────────


class TestGroupPathsByTarget:
    def test_matching_path_is_grouped(self, tmp_path):
        ursa = _make_ursa(tmp_path)
        path = _make_path("CCO")
        groups = ursa._group_paths_by_target([path], ["CCO"])
        assert any(path in v for v in groups.values())

    def test_off_target_path_discarded(self, tmp_path):
        ursa = _make_ursa(tmp_path)
        path_target = _make_path("CCO")
        path_other = _make_path("CCCO")
        groups = ursa._group_paths_by_target([path_target, path_other], ["CCO"])
        all_paths = [p for paths in groups.values() for p in paths]
        assert path_other not in all_paths

    def test_target_without_path_has_empty_list(self, tmp_path):
        ursa = _make_ursa(tmp_path)
        groups = ursa._group_paths_by_target([], ["CCO"])
        assert groups[_canonical("CCO")] == []

    def test_canonical_smiles_normalisation(self, tmp_path):
        ursa = _make_ursa(tmp_path)
        # "c1ccccc1" and "C1=CC=CC=C1" are the same molecule
        path = _make_path("c1ccccc1")
        groups = ursa._group_paths_by_target([path], ["C1=CC=CC=C1"])
        all_paths = [p for paths in groups.values() for p in paths]
        assert path in all_paths

    def test_invalid_target_smiles_skipped(self, tmp_path):
        ursa = _make_ursa(tmp_path)
        groups = ursa._group_paths_by_target([], ["NOT_VALID!!!", "CCO"])
        assert len(groups) == 1  # only CCO survives

    def test_duplicate_targets_deduplicated(self, tmp_path):
        ursa = _make_ursa(tmp_path)
        # OCC canonicalises to same as CCO
        groups = ursa._group_paths_by_target([], ["CCO", "CCO", "OCC"])
        assert len(groups) == 1

    def test_order_preserved(self, tmp_path):
        ursa = _make_ursa(tmp_path)
        targets = ["CCO", "CCCO", "CCCCO"]
        groups = ursa._group_paths_by_target([], targets)
        assert list(groups.keys()) == [_canonical(s) for s in targets]

    def test_multiple_paths_same_target_all_included(self, tmp_path):
        ursa = _make_ursa(tmp_path)
        p1 = _make_path("CCO", path_id="p1")
        p2 = _make_path("CCO", path_id="p2")
        groups = ursa._group_paths_by_target([p1, p2], ["CCO"])
        all_paths = [p for paths in groups.values() for p in paths]
        assert p1 in all_paths and p2 in all_paths


# ── SelectBestPathResult ──────────────────────────────────────────────────────


class TestSelectBestPathResult:
    def test_single_result_returned(self, tmp_path, path_1step):
        ursa = _make_ursa(tmp_path)
        pr = _make_path_result(path_1step, (3.0,))
        assert ursa._select_best_path_result([pr]) is pr

    def test_solved_preferred_over_unsolved(self, tmp_path, path_1step):
        ursa = _make_ursa(tmp_path)
        unsolved = _make_path_result(path_1step, (0.0,), is_route_solved=False)
        solved = _make_path_result(path_1step, (1.0,), is_route_solved=True)
        assert ursa._select_best_path_result([unsolved, solved]) is solved

    def test_among_solved_higher_score_wins(self, tmp_path, path_1step):
        ursa = _make_ursa(tmp_path)
        low = _make_path_result(path_1step, (1.0,), is_route_solved=True)
        high = _make_path_result(path_1step, (4.0,), is_route_solved=True)
        assert ursa._select_best_path_result([low, high]) is high

    def test_among_unsolved_better_variant_wins(self, tmp_path, path_1step):
        ursa = _make_ursa(tmp_path)
        worse = _make_path_result(path_1step, (0.0,), is_route_solved=False)
        better = _make_path_result(path_1step, (2.0,), is_route_solved=False)
        assert ursa._select_best_path_result([worse, better]) is better

    def test_all_unsolved_still_returns_a_result(self, tmp_path, path_1step):
        ursa = _make_ursa(tmp_path)
        pr1 = _make_path_result(path_1step, (0.0,), is_route_solved=False)
        pr2 = _make_path_result(path_1step, (0.5,), is_route_solved=False)
        result = ursa._select_best_path_result([pr1, pr2])
        assert result in (pr1, pr2)


# ── ScoreDataset with target_smiles ───────────────────────────────────────────


class TestScoreDatasetWithTargets:
    def test_returns_dataset_result(self, tmp_path, path_1step):
        ursa = _make_ursa(tmp_path, catalog_lines=["CC", "O"])
        result = ursa.score_dataset([path_1step], target_smiles=["CCO"])
        assert isinstance(result, DatasetResult)

    def test_total_molecules_equals_target_count(self, tmp_path, path_1step):
        ursa = _make_ursa(tmp_path, catalog_lines=["CC", "O"])
        result = ursa.score_dataset([path_1step], target_smiles=["CCO", "CCCO"])
        assert result.metrics.total_molecules == 2

    def test_off_target_paths_discarded(self, tmp_path, path_1step):
        ursa = _make_ursa(tmp_path, catalog_lines=["CC", "O"])
        off_target = _make_path("CCCO")
        result = ursa.score_dataset([path_1step, off_target], target_smiles=["CCO"])
        assert len(result.path_results) == 1

    def test_target_without_path_not_in_path_results(self, tmp_path, path_1step):
        ursa = _make_ursa(tmp_path, catalog_lines=["CC", "O"])
        # path_1step root is CCO; CCCO has no path → only 1 PathResult
        result = ursa.score_dataset([path_1step], target_smiles=["CCO", "CCCO"])
        assert len(result.path_results) == 1
        assert result.metrics.total_molecules == 2

    def test_solv_2_denominator_is_target_count(self, tmp_path, path_1step):
        ursa = _make_ursa(tmp_path, catalog_lines=["CC", "O"])
        result = ursa.score_dataset([path_1step], target_smiles=["CCO", "CCCO"])
        assert result.metrics.solv_2 == pytest.approx(result.metrics.solved_routes / 2)

    def test_both_args_raises(self, tmp_path, path_1step):
        ursa = _make_ursa(tmp_path, catalog_lines=["CC", "O"])
        with pytest.raises(ValueError, match="not both"):
            ursa.score_dataset(
                [path_1step],
                target_smiles=["CCO"],
                total_molecules=5,
            )

    def test_multiple_paths_same_target_one_result(self, tmp_path):
        ursa = _make_ursa(tmp_path, catalog_lines=["CC"])
        p1 = _make_path("CCO", path_id="p1")
        p2 = _make_path("CCO", path_id="p2")
        result = ursa.score_dataset([p1, p2], target_smiles=["CCO"])
        assert len(result.path_results) == 1

    def test_canonical_target_matching(self, tmp_path, path_1step):
        ursa = _make_ursa(tmp_path, catalog_lines=["CC", "O"])
        # OCC canonicalises to same SMILES as CCO (path_1step root)
        result = ursa.score_dataset([path_1step], target_smiles=["OCC"])
        assert len(result.path_results) == 1

    def test_empty_paths_all_targets_missing(self, tmp_path):
        ursa = _make_ursa(tmp_path)
        result = ursa.score_dataset([], target_smiles=["CCO", "CCCO"])
        assert len(result.path_results) == 0
        assert result.metrics.total_molecules == 2
        assert result.metrics.solv_2 == pytest.approx(0.0)

    def test_empty_targets_empty_result(self, tmp_path, path_1step):
        ursa = _make_ursa(tmp_path, catalog_lines=["CC", "O"])
        result = ursa.score_dataset([path_1step], target_smiles=[])
        assert len(result.path_results) == 0
        assert result.metrics.total_molecules == 0


# ── ScoreDataset backward compatibility ──────────────────────────────────────


class TestScoreDatasetBackwardCompat:
    def test_no_targets_scores_all_paths(self, tmp_path, path_1step):
        ursa = _make_ursa(tmp_path, catalog_lines=["CC", "O"])
        result = ursa.score_dataset([path_1step, path_1step])
        assert len(result.path_results) == 2

    def test_total_molecules_explicit(self, tmp_path, path_1step):
        ursa = _make_ursa(tmp_path, catalog_lines=["CC", "O"])
        result = ursa.score_dataset([path_1step], total_molecules=10)
        assert result.metrics.total_molecules == 10

    def test_total_molecules_defaults_to_len_paths(self, tmp_path, path_1step):
        ursa = _make_ursa(tmp_path, catalog_lines=["CC", "O"])
        result = ursa.score_dataset([path_1step, path_1step])
        assert result.metrics.total_molecules == 2
