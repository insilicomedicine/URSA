import pytest

from ursa.basic.node import RetrosyntheticNode
from ursa.basic.path import RetrosyntheticPath
from ursa.validation.path_consistency_checker import PathConsistencyChecker


@pytest.fixture
def checker() -> PathConsistencyChecker:
    return PathConsistencyChecker()


class TestCheck:
    def test_valid_1step_path(self, checker, path_1step):
        assert checker.check(path_1step) is True

    def test_valid_3step_linear(self, checker, path_3step_linear):
        assert checker.check(path_3step_linear) is True

    def test_valid_3step_branching(self, checker, path_3step_branching):
        assert checker.check(path_3step_branching) is True

    def test_path_with_invalid_smiles_leaf(self, checker, path_invalid_smiles):
        assert checker.check(path_invalid_smiles) is False

    def test_zero_depth_path_fails(self, checker):
        # Root is a leaf: no steps, depth=0
        root = RetrosyntheticNode(smiles="CCO")
        path = RetrosyntheticPath(path_id="empty", root=root)
        assert checker.check(path) is False

    def test_invalid_smiles_at_root_fails(self, checker):
        leaf = RetrosyntheticNode(smiles="CC")
        root = RetrosyntheticNode(smiles="NOT_VALID!!!", children=(leaf,))
        path = RetrosyntheticPath(path_id="inv_root", root=root)
        assert checker.check(path) is False

    def test_invalid_smiles_at_internal_node_fails(self, checker):
        leaf = RetrosyntheticNode(smiles="CC")
        invalid_inter = RetrosyntheticNode(smiles="NOT_VALID!!!", children=(leaf,))
        root = RetrosyntheticNode(smiles="CCO", children=(invalid_inter,))
        path = RetrosyntheticPath(path_id="inv_inter", root=root)
        assert checker.check(path) is False


class TestCheckNode:
    def test_leaf_with_valid_smiles(self, checker):
        node = RetrosyntheticNode(smiles="CCO")
        assert checker._check_node(node) is True

    def test_leaf_with_invalid_smiles(self, checker):
        node = RetrosyntheticNode(smiles="NOT_VALID!!!")
        assert checker._check_node(node) is False

    def test_internal_node_all_valid_children(self, checker):
        child = RetrosyntheticNode(smiles="CC")
        parent = RetrosyntheticNode(smiles="CCO", children=(child,))
        assert checker._check_node(parent) is True

    def test_internal_node_invalid_child(self, checker):
        child = RetrosyntheticNode(smiles="INVALID!!!")
        parent = RetrosyntheticNode(smiles="CCO", children=(child,))
        assert checker._check_node(parent) is False
