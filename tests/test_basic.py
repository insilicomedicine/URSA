import json
from pathlib import Path

import pytest
from rdkit import Chem

from ursa.basic.node import RetrosyntheticNode
from ursa.basic.path import RetrosyntheticPath

EXAMPLE_ROUTE = Path(__file__).parent / "retrocast_route_example.json"


# ── helpers ───────────────────────────────────────────────────────────────────


def _canonical(smiles: str) -> str:
    mol = Chem.MolFromSmiles(smiles)
    return Chem.MolToSmiles(mol) if mol is not None else ""


# ── fixtures: hand-built trees ────────────────────────────────────────────────


@pytest.fixture(scope="session")
def path_3step() -> RetrosyntheticPath:
    """3-step linear route.

    Tree shape:
        root
        ├── inter1
        │   └── inter2
        │       ├── leaf_boc_amine    (building block)
        │       └── leaf_iodobenzene  (building block)
        └── leaf_bromo_thiazole       (building block)
    """
    leaf_boc_amine = RetrosyntheticNode(smiles="CC(C)(C)OC(=O)N1CCn2nc(N)cc2C1")
    leaf_iodobenzene = RetrosyntheticNode(smiles="Clc1ccccc1I")
    leaf_bromo_thiazole = RetrosyntheticNode(smiles="COCc1csc(Br)n1")
    inter2 = RetrosyntheticNode(
        smiles="CC(C)(C)OC(=O)N1CCn2nc(Nc3ccccc3Cl)cc2C1",
        children=(leaf_boc_amine, leaf_iodobenzene),
    )
    inter1 = RetrosyntheticNode(
        smiles="Clc1ccccc1Nc1cc2n(n1)CCNC2",
        children=(inter2,),
    )
    root = RetrosyntheticNode(
        smiles="COCc1csc(N2CCn3nc(Nc4ccccc4Cl)cc3C2)n1",
        children=(inter1, leaf_bromo_thiazole),
    )
    return RetrosyntheticPath(path_id="X404-133-6064", root=root)


@pytest.fixture(scope="session")
def path_5step() -> RetrosyntheticPath:
    """5-step branching route: depth-4 left branch, depth-2 right branch.

    Tree shape:
        root
        ├── child1
        │   └── child1_1
        │       └── child1_1_1
        │           ├── leaf_butanoyl_cl   (building block)
        │           └── leaf_fluoroaniline (building block)
        └── child2
            ├── leaf_boronate              (building block)
            └── leaf_bromo_ether           (building block)
    """
    leaf_butanoyl_cl = RetrosyntheticNode(smiles="CCCC(=O)Cl")
    leaf_fluoroaniline = RetrosyntheticNode(smiles="NC(=O)c1ccc(F)cc1N")
    leaf_boronate = RetrosyntheticNode(smiles="CC1(C)OB(c2cc(F)c(N)c(F)c2)OC1(C)C")
    leaf_bromo_ether = RetrosyntheticNode(smiles="COc1ccc(Br)c2[nH]ncc12")
    child1_1_1 = RetrosyntheticNode(
        smiles="CCCC(=O)Nc1cc(F)ccc1C(N)=O",
        children=(leaf_butanoyl_cl, leaf_fluoroaniline),
    )
    child1_1 = RetrosyntheticNode(
        smiles="CCCc1nc(=O)c2ccc(F)cc2[nH]1",
        children=(child1_1_1,),
    )
    child1 = RetrosyntheticNode(
        smiles="CCCc1nc2cc(F)ccc2c(=O)[nH]1",
        children=(child1_1,),
    )
    child2 = RetrosyntheticNode(
        smiles="COc1ccc(-c2cc(F)c(N)c(F)c2)c2[nH]ncc12",
        children=(leaf_boronate, leaf_bromo_ether),
    )
    root = RetrosyntheticNode(
        smiles=("CCCc1nc2cc(Nc3c(F)cc(-c4ccc(OC)c5cn[nH]c45)cc3F)" "ccc2c(=O)[nH]1"),
        children=(child1, child2),
    )
    return RetrosyntheticPath(path_id="X404-133-8064", root=root)


@pytest.fixture(scope="session")
def retrocast_path() -> dict:
    """Minimal RetroCast-format dict: CCO synthesised from CC + O."""
    return {
        "target": {
            "smiles": "CCO",
            "product_of": {
                "reactants": [
                    {"smiles": "CC", "product_of": None},
                    {"smiles": "O", "product_of": None},
                ]
            },
        }
    }


# ── RetrosyntheticNode: smiles / validity ─────────────────────────────────────


@pytest.mark.parametrize(
    "smiles",
    ["CCO", "c1ccccc1", "CC(=O)O", "CN1C=NC2=C1C(=O)N(C(=O)N2C)C"],
)
def test_node_valid_smiles(smiles):
    node = RetrosyntheticNode(smiles=smiles)
    assert node.is_valid
    assert node.canonical_smiles == _canonical(smiles)


def test_node_invalid_smiles():
    node = RetrosyntheticNode(smiles="NOT_A_SMILES!!!")
    assert not node.is_valid
    assert node.canonical_smiles == ""


def test_node_canonical_normalises_equivalent_notations():
    a = RetrosyntheticNode(smiles="c1ccccc1")
    b = RetrosyntheticNode(smiles="C1=CC=CC=C1")
    assert a.canonical_smiles == b.canonical_smiles


# ── RetrosyntheticNode: structure ─────────────────────────────────────────────


def test_node_leaf_is_starting_material():
    node = RetrosyntheticNode(smiles="CCO")
    assert node.is_starting_material
    assert node.reaction_smiles == ""


def test_node_internal_is_not_starting_material():
    parent = RetrosyntheticNode(
        smiles="CCO",
        children=(RetrosyntheticNode(smiles="CC"), RetrosyntheticNode(smiles="O")),
    )
    assert not parent.is_starting_material


def test_node_reaction_smiles_uses_canonical():
    c1 = RetrosyntheticNode(smiles="CC")
    c2 = RetrosyntheticNode(smiles="O")
    parent = RetrosyntheticNode(smiles="CCO", children=(c1, c2))
    expected = f"{_canonical('CC')}.{_canonical('O')}>>{_canonical('CCO')}"
    assert parent.reaction_smiles == expected


def test_node_equality_uses_canonical():
    a = RetrosyntheticNode(smiles="c1ccccc1")
    b = RetrosyntheticNode(smiles="C1=CC=CC=C1")
    assert a == b
    assert hash(a) == hash(b)


def test_node_inequality():
    assert RetrosyntheticNode(smiles="CCO") != RetrosyntheticNode(smiles="CCCO")


# ── RetrosyntheticNode._from_mol_dict ─────────────────────────────────────────


def test_from_mol_dict_leaf_via_null_step():
    mol = {"smiles": "CCO", "product_of": None}
    node = RetrosyntheticNode._from_mol_dict(mol)
    assert node.is_starting_material
    assert node.canonical_smiles == _canonical("CCO")


def test_from_mol_dict_leaf_via_missing_step():
    mol = {"smiles": "CCO"}
    node = RetrosyntheticNode._from_mol_dict(mol)
    assert node.is_starting_material


def test_from_mol_dict_leaf_via_empty_reactants():
    mol = {"smiles": "CCO", "product_of": {"reactants": []}}
    node = RetrosyntheticNode._from_mol_dict(mol)
    assert node.is_starting_material


def test_from_mol_dict_internal_node():
    mol = {
        "smiles": "CCO",
        "product_of": {
            "reactants": [
                {"smiles": "CC", "product_of": None},
                {"smiles": "O", "product_of": None},
            ]
        },
    }
    node = RetrosyntheticNode._from_mol_dict(mol)
    assert not node.is_starting_material
    assert len(node.children) == 2
    assert node.children[0].canonical_smiles == _canonical("CC")
    assert node.children[1].canonical_smiles == _canonical("O")


# ── RetrosyntheticPath.from_dict ──────────────────────────────────────────────


def test_from_dict_path_id(retrocast_path):
    path = RetrosyntheticPath.from_dict(retrocast_path, "test-id")
    assert path.path_id == "test-id"


def test_from_dict_root_smiles(retrocast_path):
    path = RetrosyntheticPath.from_dict(retrocast_path, "p0")
    assert path.root.canonical_smiles == _canonical("CCO")


def test_from_dict_equals_manual_build(retrocast_path):
    path = RetrosyntheticPath.from_dict(retrocast_path, "p0")
    manual = RetrosyntheticPath(
        path_id="p0",
        root=RetrosyntheticNode(
            smiles="CCO",
            children=(
                RetrosyntheticNode(smiles="CC"),
                RetrosyntheticNode(smiles="O"),
            ),
        ),
    )
    assert path.root == manual.root
    assert path.num_steps == manual.num_steps
    assert path.depth == manual.depth


def test_from_dict_requires_target():
    with pytest.raises(KeyError, match="target"):
        RetrosyntheticPath.from_dict({}, path_id="x")


# ── RetrosyntheticPath: 3-step route ─────────────────────────────────────────


def test_path_3step_id(path_3step):
    assert path_3step.path_id == "X404-133-6064"


def test_path_3step_num_steps(path_3step):
    assert path_3step.num_steps == 3


def test_path_3step_depth(path_3step):
    assert path_3step.depth == 3


def test_path_3step_all_steps_dfs_order(path_3step):
    expected = [
        "COCc1csc(N2CCn3nc(Nc4ccccc4Cl)cc3C2)n1",
        "Clc1ccccc1Nc1cc2n(n1)CCNC2",
        "CC(C)(C)OC(=O)N1CCn2nc(Nc3ccccc3Cl)cc2C1",
    ]
    steps = path_3step.get_all_steps()
    assert [n.canonical_smiles for n in steps] == expected
    assert all(not n.is_starting_material for n in steps)


def test_path_3step_starting_materials_dfs_order(path_3step):
    expected = [
        "CC(C)(C)OC(=O)N1CCn2nc(N)cc2C1",
        "Clc1ccccc1I",
        "COCc1csc(Br)n1",
    ]
    sms = path_3step.get_starting_materials()
    assert [n.canonical_smiles for n in sms] == expected
    assert all(n.is_starting_material for n in sms)


def test_path_3step_reaction_smiles(path_3step):
    steps = path_3step.get_all_steps()
    for step in steps:
        assert ">>" in step.reaction_smiles
        assert step.canonical_smiles in step.reaction_smiles


# ── RetrosyntheticPath: 5-step route ─────────────────────────────────────────


def test_path_5step_id(path_5step):
    assert path_5step.path_id == "X404-133-8064"


def test_path_5step_num_steps(path_5step):
    assert path_5step.num_steps == 5


def test_path_5step_depth(path_5step):
    assert path_5step.depth == 4


def test_path_5step_all_steps_dfs_order(path_5step):
    expected = [
        "CCCc1nc2cc(Nc3c(F)cc(-c4ccc(OC)c5cn[nH]c45)cc3F)ccc2c(=O)[nH]1",
        "CCCc1nc2cc(F)ccc2c(=O)[nH]1",
        "CCCc1nc(=O)c2ccc(F)cc2[nH]1",
        "CCCC(=O)Nc1cc(F)ccc1C(N)=O",
        "COc1ccc(-c2cc(F)c(N)c(F)c2)c2[nH]ncc12",
    ]
    steps = path_5step.get_all_steps()
    assert [n.canonical_smiles for n in steps] == expected
    assert steps[0] == path_5step.root


def test_path_5step_starting_materials_dfs_order(path_5step):
    expected = [
        "CCCC(=O)Cl",
        "NC(=O)c1ccc(F)cc1N",
        "CC1(C)OB(c2cc(F)c(N)c(F)c2)OC1(C)C",
        "COc1ccc(Br)c2[nH]ncc12",
    ]
    sms = path_5step.get_starting_materials()
    assert [n.canonical_smiles for n in sms] == expected
    assert all(n.is_starting_material for n in sms)


def test_path_5step_reaction_smiles(path_5step):
    steps = path_5step.get_all_steps()
    for step in steps:
        assert ">>" in step.reaction_smiles
        assert step.canonical_smiles in step.reaction_smiles


# ── RetrosyntheticPath: RetroCast example route ───────────────────────────────


@pytest.fixture(scope="session")
def retrocast_example() -> dict:
    return json.loads(EXAMPLE_ROUTE.read_text())


@pytest.fixture(scope="session")
def path_example(retrocast_example) -> RetrosyntheticPath:
    return RetrosyntheticPath.from_dict(retrocast_example, path_id="example")


def test_example_root_smiles(retrocast_example, path_example):
    assert path_example.root.smiles == retrocast_example["target"]["smiles"]


def test_example_depth(path_example):
    assert path_example.depth == 5


def test_example_num_steps(path_example):
    assert path_example.num_steps == 7


def _json_leaf_smiles(mol: dict) -> set[str]:
    """Collect leaf SMILES by walking a RetroCast molecule's ``product_of`` tree."""
    step = mol.get("product_of")
    reactants = step.get("reactants") if step else None
    if not reactants:
        return {mol["smiles"]}
    leaves: set[str] = set()
    for r in reactants:
        leaves |= _json_leaf_smiles(r)
    return leaves


def test_example_leaves_match_json(retrocast_example, path_example):
    leaf_smiles = {n.smiles for n in path_example.get_starting_materials()}
    expected = _json_leaf_smiles(retrocast_example["target"])
    assert leaf_smiles == expected


def test_example_all_nodes_valid(path_example):
    all_nodes = path_example.get_all_steps() + path_example.get_starting_materials()
    assert all(n.is_valid for n in all_nodes)
