import json

import pytest

from tests.conftest import make_variant
from ursa.basic.building_block import BuildingBlock
from ursa.basic.result import DatasetMetrics
from ursa.basic.result import DatasetResult
from ursa.basic.result import PathResult


def _make_dataset_result(path, scores):
    vr = make_variant(path, scores)
    pr = PathResult(
        original_path=path,
        is_consistent=True,
        starting_materials=(BuildingBlock(smiles="CC", found_in_catalog=True),),
        all_bb_found=True,
        best_variant=vr,
        is_route_solved=True,
    )
    metrics = DatasetMetrics(
        total_molecules=1,
        molecules_with_route=1,
        solved_routes=1,
        passed_steps=len(scores),
        total_steps=len(scores),
        solv_2=1.0,
        mean_chemcensor_score=sum(scores) / len(scores),
    )
    return DatasetResult(path_results=(pr,), metrics=metrics)


class TestSave:
    def test_creates_two_files(self, tmp_path, path_1step):
        dr = _make_dataset_result(path_1step, (3.0,))
        metrics_path, paths_path = dr.save(tmp_path)
        assert metrics_path.exists()
        assert paths_path.exists()

    def test_default_filenames(self, tmp_path, path_1step):
        dr = _make_dataset_result(path_1step, (3.0,))
        metrics_path, paths_path = dr.save(tmp_path)
        assert metrics_path.name == "metrics.json"
        assert paths_path.name == "best_paths.json"

    def test_stem_prefix_applied(self, tmp_path, path_1step):
        dr = _make_dataset_result(path_1step, (3.0,))
        metrics_path, paths_path = dr.save(tmp_path, stem="run1")
        assert metrics_path.name == "run1_metrics.json"
        assert paths_path.name == "run1_best_paths.json"

    def test_metrics_file_is_valid_json(self, tmp_path, path_1step):
        dr = _make_dataset_result(path_1step, (3.0,))
        metrics_path, _ = dr.save(tmp_path)
        data = json.loads(metrics_path.read_text())
        assert "solv_2" in data
        assert "solved_routes" in data

    def test_paths_file_is_list(self, tmp_path, path_1step):
        dr = _make_dataset_result(path_1step, (3.0,))
        _, paths_path = dr.save(tmp_path)
        data = json.loads(paths_path.read_text())
        assert isinstance(data, list)
        assert len(data) == 1

    def test_paths_file_contains_step_scores(self, tmp_path, path_1step):
        dr = _make_dataset_result(path_1step, (3.0,))
        _, paths_path = dr.save(tmp_path)
        data = json.loads(paths_path.read_text())
        steps = data[0]["best_variant"]["steps"]
        assert len(steps) == 1
        assert steps[0]["score"] == pytest.approx(3.0)

    def test_creates_output_dir_if_missing(self, tmp_path, path_1step):
        dr = _make_dataset_result(path_1step, (3.0,))
        new_dir = tmp_path / "nested" / "output"
        dr.save(new_dir)
        assert new_dir.exists()

    def test_returns_correct_paths(self, tmp_path, path_1step):
        dr = _make_dataset_result(path_1step, (3.0,))
        metrics_path, paths_path = dr.save(tmp_path, stem="x")
        assert metrics_path == tmp_path / "x_metrics.json"
        assert paths_path == tmp_path / "x_best_paths.json"
