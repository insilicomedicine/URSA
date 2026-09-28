from pathlib import Path
from unittest.mock import patch

from ursa.data.chemcensor_db import ensure_chemcensor_db


class TestEnsureChemcensorDb:
    def test_returns_existing_target(self, tmp_path: Path) -> None:
        db_path = tmp_path / "ChemCensor-DB-U3.sqlite"
        db_path.write_bytes(b"existing")

        assert ensure_chemcensor_db(db_path) == db_path

    @patch("ursa.data.chemcensor_db.hf_hub_download")
    def test_ignores_incompatible_databases(
        self, mock_download, tmp_path: Path
    ) -> None:
        (tmp_path / "ChemCensor_DB_v1_0_0.sqlite").write_bytes(b"legacy")
        (tmp_path / "ChemCensor-DB-U2-1.0.0.sqlite").write_bytes(b"u2")
        target = tmp_path / "ChemCensor-DB-U3.sqlite"

        def download(**kwargs):
            path = tmp_path / kwargs["filename"]
            path.write_bytes(b"u3")
            return str(path)

        mock_download.side_effect = download

        assert ensure_chemcensor_db(target) == target
        assert target.read_bytes() == b"u3"
        mock_download.assert_called_once()
        assert mock_download.call_args.kwargs["filename"] == "ChemCensor-DB-U3.sqlite"
        assert (
            mock_download.call_args.kwargs["revision"]
            == "2ad529d5fd5850ce11ba6bb409e4f99a1f735cbe"
        )

    @patch("ursa.data.chemcensor_db.hf_hub_download")
    def test_downloads_sqlite_when_missing(self, mock_download, tmp_path: Path) -> None:
        target = tmp_path / "ChemCensor-DB-U3.sqlite"
        downloaded = tmp_path / "cached.sqlite"
        downloaded.write_bytes(b"downloaded")
        mock_download.return_value = str(downloaded)

        result = ensure_chemcensor_db(target)

        assert result == target
        assert target.read_bytes() == b"downloaded"
        assert not downloaded.exists()
        mock_download.assert_called_once()
