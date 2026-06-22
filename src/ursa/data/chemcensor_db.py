from __future__ import annotations

import logging
import shutil
import zipfile
from pathlib import Path

from huggingface_hub import hf_hub_download

from ..configs.data_config import DataConfig

logger = logging.getLogger(__name__)

_HF_REPO_ID = "insilicomedicine/chemcensor"
_HF_ARCHIVE = "ChemCensor-DB-U2-1.0.0.sqlite.zip"
_HF_SQLITE_NAME = "ChemCensor-DB-U2-1.0.0.sqlite"
_LEGACY_SQLITE_NAME = "ChemCensor_DB_v1_0_0.sqlite"


def ensure_chemcensor_db(db_path: Path | None = None) -> Path:
    """Return a local ChemCensor SQLite database, downloading it if needed.

    Checks the configured path first, then legacy filenames in the same
    directory. When nothing is found, downloads
    ``ChemCensor-DB-U2-1.0.0.sqlite.zip`` from HuggingFace and extracts the
    SQLite file to :attr:`~ursa.DataConfig.chemcensor_db_path`.

    :param db_path: Target database path. Defaults to
        :attr:`~ursa.DataConfig.chemcensor_db_path`.
    :return: Path to an existing SQLite database file.
    :raises RuntimeError: If the HuggingFace archive does not contain a
        ``.sqlite`` file.
    """
    target = Path(db_path) if db_path is not None else DataConfig.chemcensor_db_path
    existing = _find_existing_db(target)
    if existing is not None:
        return existing

    target.parent.mkdir(parents=True, exist_ok=True)
    logger.info("Downloading ChemCensor database from HuggingFace (%s)...", _HF_REPO_ID)
    zip_path = Path(
        hf_hub_download(
            repo_id=_HF_REPO_ID,
            repo_type="dataset",
            filename=_HF_ARCHIVE,
            revision="1eaa7a18de006147c5d1bc904b3632e1cb92420c",
            local_dir=target.parent,
        )
    )
    _extract_sqlite(zip_path, target)
    logger.info("ChemCensor database ready at %s", target)
    return target


def _find_existing_db(target: Path) -> Path | None:
    for candidate in (
        target,
        target.parent / _LEGACY_SQLITE_NAME,
        target.parent / _HF_SQLITE_NAME,
    ):
        if candidate.is_file():
            return candidate
    return None


def _extract_sqlite(zip_path: Path, target: Path) -> None:
    with zipfile.ZipFile(zip_path) as archive:
        sqlite_members = [
            name for name in archive.namelist() if name.endswith(".sqlite")
        ]
        if len(sqlite_members) != 1:
            raise RuntimeError(
                f"Expected one .sqlite file in {_HF_ARCHIVE}, "
                f"found {len(sqlite_members)}"
            )
        member = sqlite_members[0]
        with archive.open(member) as src, target.open("wb") as dst:
            shutil.copyfileobj(src, dst)
