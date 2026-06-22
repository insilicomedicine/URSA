from __future__ import annotations

import logging
import shutil
from pathlib import Path

from huggingface_hub import hf_hub_download

logger = logging.getLogger(__name__)

_HF_REPO_ID = "insilicomedicine/URSA-benchmarking-sets"
_HF_BENCHMARK_FILES = frozenset(
    {
        "URSA-expert-2026.csv",
        "URSA-drugs-clinicals-2026.csv",
    }
)


def ensure_benchmark_csv(csv_path: Path) -> Path:
    """Return a local benchmark CSV, downloading it from HuggingFace if needed.

    Only files published in ``insilicomedicine/URSA-benchmarking-sets`` are
    downloaded automatically. Custom CSV paths are returned unchanged.

    :param csv_path: Target path for the benchmark CSV.
    :return: Path to an existing CSV file.
    """
    target = Path(csv_path)
    if target.is_file():
        return target
    if target.name not in _HF_BENCHMARK_FILES:
        return target

    target.parent.mkdir(parents=True, exist_ok=True)
    logger.info(
        "Downloading benchmark set %s from HuggingFace (%s)...",
        target.name,
        _HF_REPO_ID,
    )
    downloaded = Path(
        hf_hub_download(
            repo_id=_HF_REPO_ID,
            repo_type="dataset",
            filename=target.name,
            revision="616394dac893e1cf88509aaf77fadb2a813b4226",
            local_dir=target.parent,
        )
    )
    if downloaded != target:
        shutil.move(downloaded, target)
    logger.info("Benchmark set ready at %s", target)
    return target
