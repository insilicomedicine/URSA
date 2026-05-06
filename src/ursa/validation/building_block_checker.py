from __future__ import annotations

import csv
from os import PathLike

from rdkit import Chem

from ..basic import BuildingBlock
from ..basic import RetrosyntheticPath


class BuildingBlockChecker:
    """Checks whether all leaf molecules of a path are in a catalog.

    The catalog is loaded once from a plain-text SMILES file (one SMILES
    per line) and stored as a frozen set for O(1) lookup.

    :param catalog: Set of canonical SMILES strings representing
        commercially available building blocks.
    :type catalog: frozenset[str]
    """

    def __init__(self, catalog: frozenset[str]) -> None:
        """Initialize BuildingBlockChecker with a pre-loaded catalog.

        Prefer :meth:`from_file` for constructing from a SMILES file.

        :param catalog: Frozen set of canonical building-block SMILES.
        :type catalog: frozenset[str]
        """
        self._catalog = catalog

    @classmethod
    def from_file(cls, catalog_path: str | PathLike) -> BuildingBlockChecker:
        """Construct a :class:`BuildingBlockChecker` from a SMILES file.

        Each non-empty, non-comment line in the file is treated as one
        SMILES string.

        :param catalog_path: Path to the building-block SMILES catalog file.
        :type catalog_path: str | PathLike

        :return: A new :class:`BuildingBlockChecker` instance.
        :rtype: BuildingBlockChecker
        """
        obj = cls(frozenset())
        catalog = obj._load_catalog(catalog_path)
        obj._catalog = catalog
        return obj

    def check(self, path: RetrosyntheticPath) -> tuple[BuildingBlock, ...]:
        """Check all starting-material (leaf) nodes of ``path`` against the catalog.

        Collects all nodes where
        :attr:`~ursa.RetrosyntheticNode.is_starting_material` is ``True``
        via :meth:`~ursa.RetrosyntheticPath.get_starting_materials` and looks
        up each SMILES in the catalog.

        :param path: The retrosynthetic path whose starting materials to check.
        :type path: RetrosyntheticPath

        :return: Tuple of :class:`~ursa.BuildingBlock` objects, one per
            starting-material node, with ``found_in_catalog`` set accordingly.
        :rtype: tuple[BuildingBlock, ...]
        """
        return tuple(
            BuildingBlock(
                smiles=node.smiles,
                found_in_catalog=node.canonical_smiles in self._catalog,
            )
            for node in path.get_starting_materials()
        )

    def _load_catalog(self, catalog_path: str | PathLike) -> frozenset[str]:
        """Load and return a frozen set of canonical SMILES from a catalog file.

        Accepts two formats:
        - **CSV** with a ``smiles`` column header (e.g. ``bb_id,smiles``).
        - **Plain text** — one SMILES per line; skips blanks and ``#`` comments.

        :param catalog_path: Path to the building-block catalog file.
        :type catalog_path: str | PathLike

        :return: Frozen set of canonical SMILES strings.
        :rtype: frozenset[str]
        """
        smiles_set: set[str] = set()
        with open(catalog_path, newline="") as f:
            first = f.readline()
            f.seek(0)
            if "smiles" in first.lower():
                reader = csv.DictReader(f)
                col = next(
                    c for c in (reader.fieldnames or []) if c.lower() == "smiles"
                )
                rows: list[str] = [row[col] for row in reader if row[col].strip()]
            else:
                rows = [
                    line.strip()
                    for line in f
                    if line.strip() and not line.startswith("#")
                ]
        for smi in rows:
            mol = Chem.MolFromSmiles(smi)
            if mol is not None:
                smiles_set.add(Chem.MolToSmiles(mol))
        return frozenset(smiles_set)
