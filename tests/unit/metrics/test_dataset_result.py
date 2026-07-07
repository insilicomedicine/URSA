import json

import pytest

from tests.conftest import make_path_result
from ursa.basic.result import DatasetMetrics
from ursa.basic.result import DatasetResult


def _make_dataset_result(path, scores):
    pr = make_path_result(path, scores)
    mean = sum(scores) / len(scores)
    metrics = DatasetMetrics(
        total_molecules=1,
        molecules_with_route=1,
        routes_no_synthesis=0,
        routes_solv_0=1,
        routes_solv_1=1,
        routes_solv_2=1,
        solv_0=1.0,
        solv_1=1.0,
        solv_2=1.0,
        mean_score_without_fg=mean,
        mean_score_with_fg=mean,
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
        assert "solv_0" in data
        assert "solv_1" in data
        assert "solv_2" in data
        assert "routes_solv_2" in data

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
        steps = data[0]["best_variant_solv_2"]["steps"]
        assert len(steps) == 1
        assert steps[0]["score_with_fg"] == pytest.approx(3.0)
        assert steps[0]["score_without_fg"] == pytest.approx(3.0)

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
