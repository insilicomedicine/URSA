import zipfile
from pathlib import Path
from unittest.mock import patch

import pytest

from ursa.data.chemcensor_db import ensure_chemcensor_db


def _write_zip(path: Path, sqlite_name: str, payload: bytes = b"sqlite") -> None:
    with zipfile.ZipFile(path, "w") as archive:
        archive.writestr(sqlite_name, payload)


class TestEnsureChemcensorDb:
    def test_returns_existing_target(self, tmp_path: Path) -> None:
        db_path = tmp_path / "ChemCensor-DB-U2-1.0.0.sqlite"
        db_path.write_bytes(b"existing")

        assert ensure_chemcensor_db(db_path) == db_path

    def test_returns_legacy_filename(self, tmp_path: Path) -> None:
        legacy = tmp_path / "ChemCensor_DB_v1_0_0.sqlite"
        legacy.write_bytes(b"legacy")
        target = tmp_path / "ChemCensor-DB-U2-1.0.0.sqlite"

        assert ensure_chemcensor_db(target) == legacy

    @patch("ursa.data.chemcensor_db.hf_hub_download")
    def test_downloads_and_extracts_when_missing(
        self, mock_download, tmp_path: Path
    ) -> None:
        target = tmp_path / "ChemCensor-DB-U2-1.0.0.sqlite"
        zip_path = tmp_path / "ChemCensor-DB-U2-1.0.0.sqlite.zip"
        _write_zip(zip_path, "ChemCensor-DB-U2-1.0.0.sqlite", b"downloaded")
        mock_download.return_value = str(zip_path)

        result = ensure_chemcensor_db(target)

        assert result == target
        assert target.read_bytes() == b"downloaded"
        mock_download.assert_called_once()

    @patch("ursa.data.chemcensor_db.hf_hub_download")
    def test_raises_when_archive_has_no_sqlite(
        self, mock_download, tmp_path: Path
    ) -> None:
        target = tmp_path / "ChemCensor-DB-U2-1.0.0.sqlite"
        zip_path = tmp_path / "ChemCensor-DB-U2-1.0.0.sqlite.zip"
        with zipfile.ZipFile(zip_path, "w") as archive:
            archive.writestr("readme.txt", "no db here")
        mock_download.return_value = str(zip_path)

        with pytest.raises(RuntimeError, match="Expected one .sqlite file"):
            ensure_chemcensor_db(target)
