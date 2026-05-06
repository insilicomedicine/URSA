import pytest

from ursa.basic.node import RetrosyntheticNode
from ursa.basic.path import RetrosyntheticPath
from ursa.collapsing.path_collapser import PathCollapser


@pytest.fixture
def collapser() -> PathCollapser:
    return PathCollapser()


class TestCollapse:
    def test_original_always_included(self, collapser, path_3step_linear):
        variants = collapser.collapse(path_3step_linear)
        assert path_3step_linear in variants

    def test_original_is_first(self, collapser, path_3step_linear):
        variants = collapser.collapse(path_3step_linear)
        assert variants[0] == path_3step_linear

    def test_1step_path_returns_only_original(self, collapser, path_1step):
        variants = collapser.collapse(path_1step)
        assert len(variants) == 1
        assert variants[0] == path_1step

    def test_linear_3step_generates_3_variants(self, collapser, path_3step_linear):
        # original + collapse-inter1 + collapse-inter2
        variants = collapser.collapse(path_3step_linear)
        assert len(variants) == 3

    def test_branching_3step_generates_4_variants(
        self, collapser, path_3step_branching
    ):
        # original + {inter1} + {inter2} + {inter1, inter2}
        variants = collapser.collapse(path_3step_branching)
        assert len(variants) == 4

    def test_collapsed_variants_have_fewer_steps(self, collapser, path_3step_linear):
        variants = collapser.collapse(path_3step_linear)
        original_steps = variants[0].num_steps
        for v in variants[1:]:
            assert v.num_steps < original_steps

    def test_no_duplicates(self, collapser, path_3step_linear):
        variants = collapser.collapse(path_3step_linear)
        assert len(variants) == len(set(variants))

    def test_step_limit_skips_collapsing(self):
        # Build a path with 3 steps but limit set to 2
        collapser_limited = PathCollapser(combination_step_limit=2)
        leaf1 = RetrosyntheticNode(smiles="CC")
        leaf2 = RetrosyntheticNode(smiles="O")
        inter2 = RetrosyntheticNode(smiles="CCO", children=(leaf1, leaf2))
        inter1 = RetrosyntheticNode(smiles="CCOC", children=(inter2,))
        root = RetrosyntheticNode(smiles="CCOCC", children=(inter1,))
        path = RetrosyntheticPath(path_id="big", root=root)
        variants = collapser_limited.collapse(path)
        assert variants == (path,)

    def test_variants_preserve_path_id(self, collapser, path_3step_linear):
        variants = collapser.collapse(path_3step_linear)
        for v in variants:
            assert v.path_id == path_3step_linear.path_id

    def test_full_collapse_removes_all_intermediates(
        self, collapser, path_3step_branching
    ):
        # Collapsing {inter1, inter2} should give root with 4 leaf children
        variants = collapser.collapse(path_3step_branching)
        fully_collapsed = min(variants, key=lambda v: v.num_steps)
        assert fully_collapsed.num_steps == 1
        assert len(fully_collapsed.root.children) == 4


class TestFindSkippableNodes:
    def test_root_not_skippable(self, collapser, path_3step_linear):
        skippable = collapser._find_skippable_nodes(path_3step_linear)
        assert path_3step_linear.root not in skippable

    def test_all_non_root_internals_are_skippable(self, collapser, path_3step_linear):
        skippable = collapser._find_skippable_nodes(path_3step_linear)
        steps = path_3step_linear.get_all_steps()
        expected = tuple(n for n in steps if n is not path_3step_linear.root)
        assert skippable == expected

    def test_1step_has_no_skippable(self, collapser, path_1step):
        assert collapser._find_skippable_nodes(path_1step) == ()


class TestBuildAdjacency:
    def test_parent_child_are_adjacent(self, collapser, path_3step_linear):
        adj = collapser._build_adjacency(path_3step_linear)
        steps = path_3step_linear.get_all_steps()
        # inter1 is child of root; inter2 is child of inter1
        root_smiles = steps[0].canonical_smiles
        inter1_smiles = steps[1].canonical_smiles
        inter2_smiles = steps[2].canonical_smiles
        assert inter1_smiles in adj.get(root_smiles, [])
        assert root_smiles in adj.get(inter1_smiles, [])
        assert inter2_smiles in adj.get(inter1_smiles, [])

    def test_siblings_not_adjacent(self, collapser, path_3step_branching):
        adj = collapser._build_adjacency(path_3step_branching)
        steps = path_3step_branching.get_all_steps()
        # inter1 and inter2 are siblings (both children of root)
        inter1_smiles = steps[1].canonical_smiles
        inter2_smiles = steps[2].canonical_smiles
        assert inter2_smiles not in adj.get(inter1_smiles, [])

    def test_leaves_not_in_adjacency(self, collapser, path_3step_linear):
        adj = collapser._build_adjacency(path_3step_linear)
        leaves = path_3step_linear.get_starting_materials()
        for leaf in leaves:
            assert leaf.canonical_smiles not in adj


class TestGenerateValidCombinations:
    def test_linear_no_adjacent_combination(self, collapser, path_3step_linear):
        skippable = collapser._find_skippable_nodes(path_3step_linear)
        adj = collapser._build_adjacency(path_3step_linear)
        combos = collapser._generate_valid_combinations(skippable, adj)
        # {inter1, inter2} is invalid (adjacent), so max combo size is 1
        assert all(len(c) == 1 for c in combos)

    def test_branching_allows_sibling_combination(
        self, collapser, path_3step_branching
    ):
        skippable = collapser._find_skippable_nodes(path_3step_branching)
        adj = collapser._build_adjacency(path_3step_branching)
        combos = collapser._generate_valid_combinations(skippable, adj)
        # {inter1, inter2} is valid (siblings are not adjacent)
        assert any(len(c) == 2 for c in combos)

    def test_empty_skippable_returns_empty(self, collapser, path_1step):
        skippable = collapser._find_skippable_nodes(path_1step)
        adj = collapser._build_adjacency(path_1step)
        combos = collapser._generate_valid_combinations(skippable, adj)
        assert combos == ()
