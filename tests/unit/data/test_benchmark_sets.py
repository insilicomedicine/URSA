from pathlib import Path
from unittest.mock import patch

from ursa.data.benchmark_sets import ensure_benchmark_csv


class TestEnsureBenchmarkCsv:
    def test_returns_existing_target(self, tmp_path: Path) -> None:
        csv_path = tmp_path / "URSA-expert-2026.csv"
        csv_path.write_text("Structure ID,SMILES\nid1,CC\n", encoding="utf-8")

        assert ensure_benchmark_csv(csv_path) == csv_path

    def test_returns_unknown_path_without_download(self, tmp_path: Path) -> None:
        csv_path = tmp_path / "custom-benchmark.csv"

        assert ensure_benchmark_csv(csv_path) == csv_path

    @patch("ursa.data.benchmark_sets.hf_hub_download")
    def test_downloads_known_benchmark_when_missing(
        self, mock_download, tmp_path: Path
    ) -> None:
        target = tmp_path / "URSA-expert-2026.csv"

        def _fake_download(**kwargs):
            path = Path(kwargs["local_dir"]) / kwargs["filename"]
            path.parent.mkdir(parents=True, exist_ok=True)
            path.write_text("Structure ID,SMILES\nid1,CC\n", encoding="utf-8")
            return str(path)

        mock_download.side_effect = _fake_download

        result = ensure_benchmark_csv(target)

        assert result == target
        assert target.read_text(encoding="utf-8").startswith("Structure ID,SMILES")
        mock_download.assert_called_once_with(
            repo_id="insilicomedicine/URSA-benchmarking-sets",
            repo_type="dataset",
            filename="URSA-expert-2026.csv",
            revision="616394dac893e1cf88509aaf77fadb2a813b4226",
            local_dir=target.parent,
        )
