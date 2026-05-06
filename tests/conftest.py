from pathlib import Path
from unittest.mock import MagicMock

import pytest

from ursa.basic.node import RetrosyntheticNode
from ursa.basic.path import RetrosyntheticPath
from ursa.basic.result import StepResult
from ursa.basic.result import VariantResult


# ── tree builders ─────────────────────────────────────────────────────────────


def make_leaf(smiles: str) -> RetrosyntheticNode:
    return RetrosyntheticNode(smiles=smiles)


def make_internal(smiles: str, *children: RetrosyntheticNode) -> RetrosyntheticNode:
    return RetrosyntheticNode(smiles=smiles, children=children)


# ── path fixtures ─────────────────────────────────────────────────────────────


@pytest.fixture(scope="session")
def path_1step() -> RetrosyntheticPath:
    """Single-step path: root -> [leafA, leafB]."""
    root = make_internal("CCO", make_leaf("CC"), make_leaf("O"))
    return RetrosyntheticPath(path_id="p1", root=root)


@pytest.fixture(scope="session")
def path_3step_linear() -> RetrosyntheticPath:
    """3-step linear chain: root -> inter1 -> inter2 -> [leafA, leafB].

    Skippable: inter1, inter2 (adjacent to each other).
    Valid collapse combinations: {inter1}, {inter2}.
    Expected variants: 3 (original + 2 collapsed).
    """
    leafA = make_leaf("CC")
    leafB = make_leaf("O")
    inter2 = make_internal("CCO", leafA, leafB)
    inter1 = make_internal("CCOC", inter2)
    root = make_internal("CCOCC", inter1)
    return RetrosyntheticPath(path_id="p3", root=root)


@pytest.fixture(scope="session")
def path_3step_branching() -> RetrosyntheticPath:
    """3-step branching: root -> [inter1, inter2], each with 2 leaves.

    Skippable: inter1, inter2 (siblings, not adjacent).
    Valid combinations: {inter1}, {inter2}, {inter1, inter2}.
    Expected variants: 4 (original + 3 collapsed).
    """
    inter1 = make_internal("CC", make_leaf("C"), make_leaf("C"))
    inter2 = make_internal("OO", make_leaf("O"), make_leaf("O"))
    root = make_internal("CCOO", inter1, inter2)
    return RetrosyntheticPath(path_id="p3b", root=root)


@pytest.fixture(scope="session")
def path_invalid_smiles() -> RetrosyntheticPath:
    """Path containing one node with an invalid SMILES."""
    invalid = RetrosyntheticNode(smiles="NOT_VALID!!!")
    root = make_internal("CCO", invalid, make_leaf("O"))
    return RetrosyntheticPath(path_id="pinv", root=root)


# ── building-block catalog fixtures ──────────────────────────────────────────


@pytest.fixture
def catalog_file(tmp_path: Path) -> Path:
    """Temp file with 3 valid SMILES and 1 comment line."""
    p = tmp_path / "catalog.smi"
    p.write_text("# comment\nCC\nO\nc1ccccc1\n\n")
    return p


@pytest.fixture
def empty_catalog_file(tmp_path: Path) -> Path:
    p = tmp_path / "empty.smi"
    p.write_text("# only comment\n")
    return p


# ── mock scorer ───────────────────────────────────────────────────────────────


def make_scorer(*scores: float):
    """Return a mock ReactionScorer that yields ``scores`` in order."""
    mock = MagicMock()
    mock.score.side_effect = list(scores)
    return mock


def make_const_scorer(value: float):
    """Return a mock ReactionScorer that always returns ``value``."""
    mock = MagicMock()
    mock.score.return_value = value
    return mock


# ── VariantResult builder ─────────────────────────────────────────────────────


def make_variant(
    path: RetrosyntheticPath,
    scores: tuple[float, ...],
    pass_threshold: float = 1.0,
) -> VariantResult:
    nodes = path.get_all_steps()
    step_results = tuple(
        StepResult(node=n, score=s, passed=s >= pass_threshold)
        for n, s in zip(nodes, scores)
    )
    total = len(step_results)
    passed = sum(1 for sr in step_results if sr.passed)
    pf = passed / total if total else 0.0
    avg = sum(sr.score for sr in step_results) / total if total else 0.0
    return VariantResult(
        path=path,
        step_results=step_results,
        percent_found=pf,
        chemcensor_per_route=avg,
        all_steps_passed=all(sr.passed for sr in step_results),
    )
