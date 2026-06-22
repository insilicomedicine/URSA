import csv
from pathlib import Path
from unittest.mock import patch

import pytest

from ursa.datasets import BenchmarkDataset
from ursa.datasets import TargetEntry


# ── helpers ───────────────────────────────────────────────────────────────────


def _write_csv(
    path: Path,
    rows: list[dict],
    fieldnames: list[str] | None = None,
) -> None:
    if fieldnames is None:
        fieldnames = list(rows[0].keys()) if rows else []
    with open(path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


# ── TestLoad ──────────────────────────────────────────────────────────────────


class TestLoad:
    def test_returns_tuple_of_target_entries(self, tmp_path):
        csv_path = tmp_path / "ds.csv"
        _write_csv(csv_path, [{"Structure ID": "m1", "SMILES": "CCO"}])
        ds = BenchmarkDataset.from_csv(csv_path)
        result = ds.load()
        assert isinstance(result, tuple)
        assert all(isinstance(e, TargetEntry) for e in result)

    def test_entry_count_matches_rows(self, tmp_path):
        csv_path = tmp_path / "ds.csv"
        _write_csv(
            csv_path,
            [
                {"Structure ID": "m1", "SMILES": "CCO"},
                {"Structure ID": "m2", "SMILES": "CCCO"},
                {"Structure ID": "m3", "SMILES": "c1ccccc1"},
            ],
        )
        ds = BenchmarkDataset.from_csv(csv_path)
        assert len(ds.load()) == 3

    def test_entry_id_and_smiles(self, tmp_path):
        csv_path = tmp_path / "ds.csv"
        _write_csv(csv_path, [{"Structure ID": "mol-1", "SMILES": "CC"}])
        entry = BenchmarkDataset.from_csv(csv_path).load()[0]
        assert entry.id == "mol-1"
        assert entry.smiles == "CC"

    def test_empty_csv_returns_empty_tuple(self, tmp_path):
        csv_path = tmp_path / "empty.csv"
        _write_csv(csv_path, [], fieldnames=["Structure ID", "SMILES"])
        ds = BenchmarkDataset.from_csv(csv_path)
        assert ds.load() == ()

    def test_custom_column_names(self, tmp_path):
        csv_path = tmp_path / "custom.csv"
        _write_csv(
            csv_path, [{"mol_id": "x1", "smi": "CC"}], fieldnames=["mol_id", "smi"]
        )
        ds = BenchmarkDataset.from_csv(csv_path, id_col="mol_id", smiles_col="smi")
        entry = ds.load()[0]
        assert entry.id == "x1"
        assert entry.smiles == "CC"

    def test_target_entry_is_frozen(self, tmp_path):
        csv_path = tmp_path / "ds.csv"
        _write_csv(csv_path, [{"Structure ID": "m1", "SMILES": "CCO"}])
        entry = BenchmarkDataset.from_csv(csv_path).load()[0]
        with pytest.raises((AttributeError, TypeError)):
            entry.id = "changed"  # type: ignore[misc]


# ── TestFromCsv ───────────────────────────────────────────────────────────────


class TestFromCsv:
    def test_default_column_names(self, tmp_path):
        csv_path = tmp_path / "ds.csv"
        _write_csv(csv_path, [{"Structure ID": "id1", "SMILES": "O"}])
        ds = BenchmarkDataset.from_csv(csv_path)
        assert ds.load()[0].id == "id1"

    def test_returns_benchmark_dataset(self, tmp_path):
        csv_path = tmp_path / "ds.csv"
        _write_csv(csv_path, [], fieldnames=["Structure ID", "SMILES"])
        assert isinstance(BenchmarkDataset.from_csv(csv_path), BenchmarkDataset)

    def test_path_as_string(self, tmp_path):
        csv_path = tmp_path / "ds.csv"
        _write_csv(csv_path, [{"Structure ID": "s1", "SMILES": "CC"}])
        ds = BenchmarkDataset.from_csv(str(csv_path))
        assert len(ds.load()) == 1


# ── TestTargetSmiles ──────────────────────────────────────────────────────────


class TestTargetSmiles:
    def test_returns_tuple_of_strings(self, tmp_path):
        csv_path = tmp_path / "ds.csv"
        _write_csv(
            csv_path,
            [
                {"Structure ID": "m1", "SMILES": "CCO"},
                {"Structure ID": "m2", "SMILES": "O"},
            ],
        )
        ds = BenchmarkDataset.from_csv(csv_path)
        smiles = ds.target_smiles
        assert isinstance(smiles, tuple)
        assert all(isinstance(s, str) for s in smiles)

    def test_smiles_order_preserved(self, tmp_path):
        csv_path = tmp_path / "ds.csv"
        _write_csv(
            csv_path,
            [
                {"Structure ID": "m1", "SMILES": "CCO"},
                {"Structure ID": "m2", "SMILES": "CCCO"},
                {"Structure ID": "m3", "SMILES": "c1ccccc1"},
            ],
        )
        ds = BenchmarkDataset.from_csv(csv_path)
        assert ds.target_smiles == ("CCO", "CCCO", "c1ccccc1")

    def test_empty_dataset_returns_empty_tuple(self, tmp_path):
        csv_path = tmp_path / "ds.csv"
        _write_csv(csv_path, [], fieldnames=["Structure ID", "SMILES"])
        assert BenchmarkDataset.from_csv(csv_path).target_smiles == ()


def _write_benchmark_csv(path: Path, rows: int) -> None:
    with open(path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=["Structure ID", "SMILES"])
        writer.writeheader()
        for i in range(rows):
            writer.writerow({"Structure ID": f"id{i}", "SMILES": "CC"})


# ── TestBundledPresets ────────────────────────────────────────────────────────


class TestBundledPresets:
    @pytest.mark.parametrize(
        "attr, expected_len",
        [
            ("EXPERT_2026", 100),
            ("DRUGS_CLINICALS_2026", 100),
        ],
    )
    def test_preset_is_benchmark_dataset(self, attr, expected_len):
        ds = getattr(BenchmarkDataset, attr)
        assert isinstance(ds, BenchmarkDataset)

    @pytest.mark.parametrize(
        "attr, expected_len",
        [
            ("EXPERT_2026", 100),
            ("DRUGS_CLINICALS_2026", 100),
        ],
    )
    @patch("ursa.datasets.benchmark_dataset.ensure_benchmark_csv")
    def test_preset_loads_correct_row_count(
        self, mock_ensure, attr, expected_len, tmp_path: Path
    ):
        ds = getattr(BenchmarkDataset, attr)
        csv_path = tmp_path / ds._path.name
        _write_benchmark_csv(csv_path, expected_len)
        mock_ensure.return_value = csv_path
        assert len(ds.load()) == expected_len

    @pytest.mark.parametrize("attr", ["EXPERT_2026", "DRUGS_CLINICALS_2026"])
    @patch("ursa.datasets.benchmark_dataset.ensure_benchmark_csv")
    def test_preset_entries_have_non_empty_ids_and_smiles(
        self, mock_ensure, attr, tmp_path: Path
    ):
        ds = getattr(BenchmarkDataset, attr)
        csv_path = tmp_path / ds._path.name
        _write_benchmark_csv(csv_path, 3)
        mock_ensure.return_value = csv_path
        for entry in ds.load():
            assert entry.id
            assert entry.smiles

    @patch("ursa.datasets.benchmark_dataset.ensure_benchmark_csv")
    def test_preset_target_smiles_length(self, mock_ensure, tmp_path: Path):
        ds = BenchmarkDataset.EXPERT_2026
        csv_path = tmp_path / ds._path.name
        _write_benchmark_csv(csv_path, 100)
        mock_ensure.return_value = csv_path
        assert len(ds.target_smiles) == 100
