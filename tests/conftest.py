from pathlib import Path
from unittest.mock import MagicMock

import pytest

from ursa.basic.building_block import BuildingBlock
from ursa.basic.node import RetrosyntheticNode
from ursa.basic.path import RetrosyntheticPath
from ursa.basic.result import PathResult
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


class _Score:
    """Minimal stand-in for ``chemcensor.ScoreResult`` used in tests."""

    def __init__(self, without_functional_groups: float, with_functional_groups: float):
        self.without_functional_groups = without_functional_groups
        self.with_functional_groups = with_functional_groups


def _to_score(value: float | tuple[float, float]) -> _Score:
    """Coerce a float (uniform) or ``(without_fg, with_fg)`` tuple to ``_Score``."""
    if isinstance(value, tuple):
        without_fg, with_fg = value
        return _Score(without_fg, with_fg)
    return _Score(value, value)


def make_scorer(*scores: float | tuple[float, float]):
    """Return a mock ReactionScorer whose ``evaluate`` yields ``scores`` in order.

    Each score is either a float (same with/without FG) or a
    ``(without_fg, with_fg)`` tuple to exercise the Solv-1/Solv-2 split.
    """
    mock = MagicMock()
    mock.evaluate.side_effect = [_to_score(s) for s in scores]
    return mock


def make_const_scorer(value: float | tuple[float, float]):
    """Return a mock ReactionScorer whose ``evaluate`` always returns ``value``."""
    mock = MagicMock()
    mock.evaluate.return_value = _to_score(value)
    return mock


# ── VariantResult builder ─────────────────────────────────────────────────────


def make_variant(
    path: RetrosyntheticPath,
    scores: tuple[float, ...],
    scores_with_fg: tuple[float, ...] | None = None,
) -> VariantResult:
    """Build a VariantResult; ``scores`` is the without-FG dimension.

    When ``scores_with_fg`` is omitted the with-FG scores mirror ``scores``.
    """
    nodes = path.get_all_steps()
    with_fg = scores_with_fg if scores_with_fg is not None else scores
    step_results = tuple(
        StepResult(node=n, score_without_fg=wo, score_with_fg=w)
        for n, wo, w in zip(nodes, scores, with_fg)
    )
    total = len(step_results)
    has_steps = total > 0
    pass_1 = has_steps and all(sr.passes_solv_1 for sr in step_results)
    pass_2 = has_steps and all(sr.passes_solv_2 for sr in step_results)
    mean_wo = (
        sum(sr.score_without_fg for sr in step_results) / total if has_steps else 0.0
    )
    mean_w = sum(sr.score_with_fg for sr in step_results) / total if has_steps else 0.0
    return VariantResult(
        path=path,
        step_results=step_results,
        all_steps_pass_solv_1=pass_1,
        all_steps_pass_solv_2=pass_2,
        mean_score_without_fg=mean_wo,
        mean_score_with_fg=mean_w,
    )


# ── PathResult builder ────────────────────────────────────────────────────────


def make_path_result(
    path: RetrosyntheticPath,
    scores: tuple[float, ...],
    *,
    scores_with_fg: tuple[float, ...] | None = None,
    bb_found: bool = True,
    is_consistent: bool = True,
    is_no_synthesis: bool | None = None,
) -> PathResult:
    """Build a PathResult using a single variant for both Solv levels.

    ``scores`` is the without-FG dimension; ``scores_with_fg`` defaults to
    mirror it. Pass flags are derived exactly as :meth:`Ursa.score` does.
    """
    variant = make_variant(path, scores, scores_with_fg)
    if is_no_synthesis is None:
        is_no_synthesis = path.num_steps == 0
    passes_solv_0 = not is_no_synthesis and is_consistent and bb_found
    passes_solv_1 = passes_solv_0 and variant.all_steps_pass_solv_1
    passes_solv_2 = passes_solv_0 and variant.all_steps_pass_solv_2
    return PathResult(
        original_path=path,
        is_consistent=is_consistent,
        starting_materials=(BuildingBlock(smiles="CC", found_in_catalog=bb_found),),
        all_bb_found=bb_found,
        is_no_synthesis=is_no_synthesis,
        passes_solv_0=passes_solv_0,
        passes_solv_1=passes_solv_1,
        passes_solv_2=passes_solv_2,
        best_variant_solv_1=variant,
        best_variant_solv_2=variant,
    )
