from pathlib import Path

_DATA_DIR = Path(__file__).parent.parent.parent.parent / "data"


class DataConfig:
    """Default filesystem paths for URSA data assets.

    All paths point to the ``data/`` directory at the project root.
    Pass them explicitly to override (e.g. in tests or custom setups).
    """

    bb_catalog_path: Path = _DATA_DIR / "building_blocks" / "URSA_BBs_v1_0_0.csv"
    chemcensor_db_path: Path = (
        _DATA_DIR / "chemcensor_db" / "ChemCensor-DB-U2-1.0.0.sqlite"
    )
    benchmark_sets_dir: Path = _DATA_DIR / "URSA_benchmarking_sets"
