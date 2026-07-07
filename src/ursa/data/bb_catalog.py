from __future__ import annotations

import logging
import shutil
from pathlib import Path

from huggingface_hub import hf_hub_download

from ..configs.data_config import DataConfig

logger = logging.getLogger(__name__)

_HF_REPO_ID = "insilicomedicine/URSA-BBs"
_HF_CATALOG_NAME = "URSA_BBs_v1_0_0.csv"


def ensure_bb_catalog(catalog_path: Path | None = None) -> Path:
    """Return a local building-block catalog, downloading it if needed.

    When the configured CSV is missing, downloads ``URSA_BBs_v1_0_0.csv`` from
    HuggingFace into the parent directory of
    :attr:`~ursa.DataConfig.bb_catalog_path`.

    :param catalog_path: Target catalog path. Defaults to
        :attr:`~ursa.DataConfig.bb_catalog_path`.
    :return: Path to an existing catalog CSV file.
    """
    target = (
        Path(catalog_path) if catalog_path is not None else DataConfig.bb_catalog_path
    )
    existing = _find_existing_catalog(target)
    if existing is not None:
        return existing

    target.parent.mkdir(parents=True, exist_ok=True)
    logger.info(
        "Downloading building-block catalog from HuggingFace (%s)...", _HF_REPO_ID
    )
    downloaded = Path(
        hf_hub_download(
            repo_id=_HF_REPO_ID,
            repo_type="dataset",
            filename=_HF_CATALOG_NAME,
            revision="9ac1948b272fafbaaae99916a02eb5298907138b",
            local_dir=target.parent,
        )
    )
    if downloaded != target:
        shutil.move(downloaded, target)
    logger.info("Building-block catalog ready at %s", target)
    return target


def _find_existing_catalog(target: Path) -> Path | None:
    if target.is_file():
        return target
    return None
