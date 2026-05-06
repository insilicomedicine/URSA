from rdkit import Chem

from ursa.basic.building_block import BuildingBlock
from ursa.validation.building_block_checker import BuildingBlockChecker


def canonical(smiles: str) -> str:
    mol = Chem.MolFromSmiles(smiles)
    return Chem.MolToSmiles(mol) if mol else ""


class TestLoadCatalog:
    def test_loads_valid_smiles(self, catalog_file):
        checker = BuildingBlockChecker.from_file(catalog_file)
        assert canonical("CC") in checker._catalog
        assert canonical("O") in checker._catalog
        assert canonical("c1ccccc1") in checker._catalog

    def test_skips_comment_lines(self, catalog_file):
        checker = BuildingBlockChecker.from_file(catalog_file)
        assert "# comment" not in checker._catalog

    def test_skips_blank_lines(self, catalog_file):
        checker = BuildingBlockChecker.from_file(catalog_file)
        assert "" not in checker._catalog

    def test_empty_catalog(self, empty_catalog_file):
        checker = BuildingBlockChecker.from_file(empty_catalog_file)
        assert len(checker._catalog) == 0

    def test_canonicalises_smiles(self, tmp_path):
        p = tmp_path / "cat.smi"
        p.write_text("C1=CC=CC=C1\n")
        checker = BuildingBlockChecker.from_file(p)
        assert canonical("c1ccccc1") in checker._catalog

    def test_skips_invalid_smiles(self, tmp_path):
        p = tmp_path / "cat.smi"
        p.write_text("CC\nNOT_VALID!!!\n")
        checker = BuildingBlockChecker.from_file(p)
        assert len(checker._catalog) == 1


class TestCheck:
    def test_all_found(self, catalog_file, path_1step):
        # path_1step leaves: CC, O — both in catalog
        checker = BuildingBlockChecker.from_file(catalog_file)
        result = checker.check(path_1step)
        assert len(result) == 2
        assert all(isinstance(bb, BuildingBlock) for bb in result)
        assert all(bb.found_in_catalog for bb in result)

    def test_none_found(self, empty_catalog_file, path_1step):
        checker = BuildingBlockChecker.from_file(empty_catalog_file)
        result = checker.check(path_1step)
        assert all(not bb.found_in_catalog for bb in result)

    def test_partial_found(self, tmp_path, path_1step):
        # path_1step leaves: CC, O — catalog only has CC
        p = tmp_path / "cat.smi"
        p.write_text("CC\n")
        checker = BuildingBlockChecker.from_file(p)
        result = checker.check(path_1step)
        found = [bb.found_in_catalog for bb in result]
        assert True in found
        assert False in found

    def test_preserves_original_smiles(self, catalog_file, path_1step):
        checker = BuildingBlockChecker.from_file(catalog_file)
        result = checker.check(path_1step)
        leaves = path_1step.get_starting_materials()
        assert tuple(bb.smiles for bb in result) == tuple(n.smiles for n in leaves)

    def test_returns_tuple(self, catalog_file, path_1step):
        checker = BuildingBlockChecker.from_file(catalog_file)
        result = checker.check(path_1step)
        assert isinstance(result, tuple)

    def test_count_matches_leaves(self, catalog_file, path_3step_linear):
        checker = BuildingBlockChecker.from_file(catalog_file)
        result = checker.check(path_3step_linear)
        assert len(result) == len(path_3step_linear.get_starting_materials())
