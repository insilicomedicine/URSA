from pathlib import Path
from unittest.mock import patch

from ursa.data.bb_catalog import ensure_bb_catalog


class TestEnsureBbCatalog:
    def test_returns_existing_target(self, tmp_path: Path) -> None:
        catalog_path = tmp_path / "URSA_BBs_v1_0_0.csv"
        catalog_path.write_text("bb_id,smiles\nURSA-BB-000000,CC\n", encoding="utf-8")

        assert ensure_bb_catalog(catalog_path) == catalog_path

    @patch("ursa.data.bb_catalog.hf_hub_download")
    def test_downloads_when_missing(self, mock_download, tmp_path: Path) -> None:
        target = tmp_path / "URSA_BBs_v1_0_0.csv"

        def _fake_download(**kwargs):
            path = Path(kwargs["local_dir"]) / kwargs["filename"]
            path.parent.mkdir(parents=True, exist_ok=True)
            path.write_text("bb_id,smiles\nURSA-BB-000000,CC\n", encoding="utf-8")
            return str(path)

        mock_download.side_effect = _fake_download

        result = ensure_bb_catalog(target)

        assert result == target
        assert target.read_text(encoding="utf-8").startswith("bb_id,smiles")
        mock_download.assert_called_once_with(
            repo_id="insilicomedicine/URSA-BBs",
            repo_type="dataset",
            filename="URSA_BBs_v1_0_0.csv",
            revision="9ac1948b272fafbaaae99916a02eb5298907138b",
            local_dir=target.parent,
        )

    @patch("ursa.data.bb_catalog.hf_hub_download")
    def test_moves_download_when_hf_returns_different_path(
        self, mock_download, tmp_path: Path
    ) -> None:
        target = tmp_path / "URSA_BBs_v1_0_0.csv"
        staging = tmp_path / ".cache" / "URSA_BBs_v1_0_0.csv"
        staging.parent.mkdir(parents=True)
        staging.write_text("bb_id,smiles\nURSA-BB-000001,CO\n", encoding="utf-8")
        mock_download.return_value = str(staging)

        result = ensure_bb_catalog(target)

        assert result == target
        assert target.read_text(encoding="utf-8").startswith("bb_id,smiles")
        assert not staging.exists()
