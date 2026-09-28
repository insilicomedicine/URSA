from __future__ import annotations

import logging
import shutil
from pathlib import Path

from huggingface_hub import hf_hub_download

from ..configs.data_config import DataConfig

logger = logging.getLogger(__name__)

_HF_REPO_ID = "insilicomedicine/chemcensor"
_HF_SQLITE_NAME = "ChemCensor-DB-U3.sqlite"


def ensure_chemcensor_db(db_path: Path | None = None) -> Path:
    """Return the ChemCensor U3 database, downloading it if needed.

    ChemCensor 1.4 requires ``ChemCensor-DB-U3.sqlite``. Older databases,
    including ``ChemCensor-DB-U2-1.0.0.sqlite`` and
    ``ChemCensor_DB_v1_0_0.sqlite``, are not accepted. When the configured
    file is missing, this downloads it from Hugging Face into
    :attr:`~ursa.DataConfig.chemcensor_db_path`.

    :param db_path: Target database path. Defaults to
        :attr:`~ursa.DataConfig.chemcensor_db_path`.
    :return: Path to the existing SQLite database file.
    """
    target = Path(db_path) if db_path is not None else DataConfig.chemcensor_db_path
    if target.is_file():
        return target

    target.parent.mkdir(parents=True, exist_ok=True)
    logger.info("Downloading ChemCensor database from HuggingFace (%s)...", _HF_REPO_ID)
    downloaded = Path(
        hf_hub_download(
            repo_id=_HF_REPO_ID,
            repo_type="dataset",
            filename=_HF_SQLITE_NAME,
            revision="2ad529d5fd5850ce11ba6bb409e4f99a1f735cbe",
            local_dir=target.parent,
        )
    )
    if downloaded != target:
        shutil.move(downloaded, target)
    logger.info("ChemCensor database ready at %s", target)
    return target
