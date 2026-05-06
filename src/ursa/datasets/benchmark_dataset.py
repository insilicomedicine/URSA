from __future__ import annotations

import csv
from dataclasses import dataclass
from pathlib import Path
from typing import ClassVar

_DATASETS_DIR = (
    Path(__file__).parent.parent.parent.parent / "data" / "URSA_benchmarking_sets"
)


@dataclass(frozen=True)
class TargetEntry:
    """A single target molecule in a benchmark dataset.

    :param id: Dataset-specific identifier (e.g. ``"X404-133-8064"``).
    :type id: str
    :param smiles: SMILES string of the target molecule.
    :type smiles: str
    """

    id: str
    smiles: str


class BenchmarkDataset:
    """Loader for target-molecule benchmark datasets.

    Supports two usage patterns:

    **Predefined bundled datasets** — select via class attributes::

        dataset = BenchmarkDataset.EXPERT_2026
        ursa.score_dataset(paths, target_smiles=dataset.target_smiles)

    **Custom CSV** — provide a path and column names::

        dataset = BenchmarkDataset.from_csv(
            "my_targets.csv", id_col="mol_id", smiles_col="smiles"
        )

    :param path: Path to a CSV file containing target molecules.
    :type path: str | Path
    :param id_col: Name of the column containing molecule identifiers.
    :type id_col: str
    :param smiles_col: Name of the column containing SMILES strings.
    :type smiles_col: str
    """

    EXPERT_2026: ClassVar[BenchmarkDataset]
    DRUGS_CLINICALS_2026: ClassVar[BenchmarkDataset]
    USPTO_190: ClassVar[BenchmarkDataset]

    def __init__(
        self,
        path: str | Path,
        *,
        id_col: str,
        smiles_col: str,
    ) -> None:
        self._path = Path(path)
        self._id_col = id_col
        self._smiles_col = smiles_col

    @classmethod
    def from_csv(
        cls,
        path: str | Path,
        *,
        id_col: str = "Structure ID",
        smiles_col: str = "SMILES",
    ) -> BenchmarkDataset:
        """Construct a :class:`BenchmarkDataset` from a CSV file.

        :param path: Path to the CSV file.
        :type path: str | Path
        :param id_col: Column name for molecule identifiers.
            Defaults to ``"Structure ID"``.
        :type id_col: str
        :param smiles_col: Column name for SMILES strings.
            Defaults to ``"SMILES"``.
        :type smiles_col: str

        :return: A new :class:`BenchmarkDataset` instance.
        :rtype: BenchmarkDataset
        """
        return cls(path, id_col=id_col, smiles_col=smiles_col)

    def load(self) -> tuple[TargetEntry, ...]:
        """Load all entries from the CSV file.

        :return: Tuple of :class:`TargetEntry` objects, one per data row.
        :rtype: tuple[TargetEntry, ...]
        """
        entries: list[TargetEntry] = []
        with open(self._path, newline="") as f:
            reader = csv.DictReader(f)
            for row in reader:
                entries.append(
                    TargetEntry(id=row[self._id_col], smiles=row[self._smiles_col])
                )
        return tuple(entries)

    @property
    def target_smiles(self) -> tuple[str, ...]:
        """SMILES strings of all targets, suitable for
        :meth:`~ursa.Ursa.score_dataset`.

        :rtype: tuple[str, ...]
        """
        return tuple(e.smiles for e in self.load())

    def __repr__(self) -> str:
        return f"BenchmarkDataset({self._path.name!r})"


# ── Bundled presets ───────────────────────────────────────────────────────────

BenchmarkDataset.EXPERT_2026 = BenchmarkDataset(
    _DATASETS_DIR / "URSA-expert-2026.csv",
    id_col="Structure ID",
    smiles_col="SMILES",
)

BenchmarkDataset.DRUGS_CLINICALS_2026 = BenchmarkDataset(
    _DATASETS_DIR / "URSA-drugs-clinicals-2026.csv",
    id_col="Structure ID",
    smiles_col="SMILES",
)

BenchmarkDataset.USPTO_190 = BenchmarkDataset(
    _DATASETS_DIR / "Uspto-190.csv",
    id_col="Structure ID",
    smiles_col="SMILES",
)
