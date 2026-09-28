import zipfile
from pathlib import Path

import pytest
from lxml import etree

from ursa.basic.node import RetrosyntheticNode
from ursa.basic.path import RetrosyntheticPath
from ursa.basic.result import DatasetMetrics
from ursa.basic.result import DatasetResult
from ursa.basic.result import PathResult
from ursa.basic.result import StepResult
from ursa.basic.result import VariantResult
from ursa.drawing import schemes
from ursa.drawing import write_schemes_zip


@pytest.fixture(autouse=True)
def _skip_chemcensor_processing(monkeypatch):
    """Keep archive tests fast; ChemCensor processing is checked separately."""
    monkeypatch.setattr(
        schemes,
        "_prepare_reactions_for_drawing",
        lambda reactions: reactions,
    )


def _path_result(path_id: str, rxn: str | None, score: float = 0.5) -> PathResult:
    """Build a ``PathResult`` from ``reactants>>product``; ``None`` → empty path."""
    if rxn is None:
        root = RetrosyntheticNode(smiles="c1ccccc1")
        path = RetrosyntheticPath(path_id=path_id, root=root)
        steps: tuple[StepResult, ...] = ()
    else:
        reactants_str, product_str = rxn.split(">>")
        leaves = tuple(
            RetrosyntheticNode(smiles=s) for s in reactants_str.split(".") if s
        )
        root = RetrosyntheticNode(smiles=product_str, children=leaves)
        path = RetrosyntheticPath(path_id=path_id, root=root)
        steps = (StepResult(node=root, score_without_fg=score, score_with_fg=score),)

    has_steps = bool(steps)
    bv = VariantResult(
        path=path,
        step_results=steps,
        all_steps_pass_solv_1=has_steps,
        all_steps_pass_solv_2=has_steps,
        mean_score_without_fg=score if has_steps else 0.0,
        mean_score_with_fg=score if has_steps else 0.0,
    )
    return PathResult(
        original_path=path,
        is_consistent=True,
        starting_materials=(),
        all_bb_found=True,
        is_no_synthesis=not has_steps,
        passes_solv_0=has_steps,
        passes_solv_1=has_steps,
        passes_solv_2=has_steps,
        best_variant_solv_1=bv,
        best_variant_solv_2=bv,
    )


def _dataset(path_results) -> DatasetResult:
    metrics = DatasetMetrics(
        total_molecules=len(path_results),
        molecules_with_route=len(path_results),
        routes_no_synthesis=0,
        routes_solv_0=0,
        routes_solv_1=0,
        routes_solv_2=0,
        solv_0=0.0,
        solv_1=0.0,
        solv_2=0.0,
        mean_score_without_fg=0.0,
        mean_score_with_fg=0.0,
    )
    return DatasetResult(path_results=tuple(path_results), metrics=metrics)


def test_write_schemes_zip_renders_stepped_paths(tmp_path: Path):
    result = _dataset([_path_result("mol1__r1", "CCO.CC(=O)O>>CCOC(C)=O", score=0.87)])

    zip_path = write_schemes_zip(result, tmp_path, stem="model")

    assert zip_path == tmp_path / "model_schemes.zip"
    with zipfile.ZipFile(zip_path) as archive:
        names = archive.namelist()
        assert names == ["mol1__r1.cdxml"]
        content = archive.read(names[0])
        # Well-formed CDXML carrying the per-step score label.
        etree.fromstring(content)
        assert b"<CDXML" in content
        assert b"Score: 0.87" in content


def test_write_schemes_zip_skips_empty_paths(tmp_path: Path):
    result = _dataset(
        [
            _path_result("withsteps__r1", "CCO.CC(=O)O>>CCOC(C)=O"),
            _path_result("empty__r1", None),
        ]
    )

    zip_path = write_schemes_zip(result, tmp_path)

    assert zip_path == tmp_path / "schemes.zip"
    with zipfile.ZipFile(zip_path) as archive:
        # Only the path that has reaction steps is rendered.
        assert archive.namelist() == ["withsteps__r1.cdxml"]


def test_write_schemes_zip_returns_none_when_nothing_rendered(tmp_path: Path):
    result = _dataset([_path_result("empty__r1", None)])

    zip_path = write_schemes_zip(result, tmp_path)

    assert zip_path is None
    assert not (tmp_path / "schemes.zip").exists()


def test_write_schemes_zip_sanitizes_and_deduplicates_filenames(tmp_path: Path):
    rxn = "CCO.CC(=O)O>>CCOC(C)=O"
    result = _dataset(
        [
            _path_result("a/b:c*?__r1", rxn),
            _path_result("a/b:c*?__r1", rxn),  # collides after sanitization
        ]
    )

    zip_path = write_schemes_zip(result, tmp_path)

    assert zip_path == tmp_path / "schemes.zip"
    with zipfile.ZipFile(zip_path) as archive:
        names = archive.namelist()
    assert names == ["a_b_c___r1.cdxml", "a_b_c___r1_2.cdxml"]


def test_write_schemes_zip_uses_target_id_for_filename(tmp_path: Path):
    result = _dataset([_path_result("CCOC(C)=O__r1", "CCO.CC(=O)O>>CCOC(C)=O")])

    zip_path = write_schemes_zip(
        result,
        tmp_path,
        target_id_by_smiles={"CCOC(C)=O": "EXPERT-324"},
    )

    assert zip_path is not None
    with zipfile.ZipFile(zip_path) as archive:
        assert archive.namelist() == ["EXPERT-324.cdxml"]


def test_write_schemes_zip_draws_processed_reactions(monkeypatch, tmp_path: Path):
    result = _dataset([_path_result("mol1__r1", "CCO.O.[Na+]>>CC=O")])
    received: list[list[str]] = []

    def prepare(reactions):
        received.append(reactions)
        return ["CCO>>CC=O"]

    monkeypatch.setattr(schemes, "_prepare_reactions_for_drawing", prepare)

    zip_path = write_schemes_zip(result, tmp_path)

    assert zip_path is not None
    assert received == [["CCO.O.[Na+]>>CC=O"]]
    with zipfile.ZipFile(zip_path) as archive:
        content = archive.read("mol1__r1.cdxml")
    assert b'Element="11"' not in content
