from __future__ import annotations

import csv
from os import PathLike

from rdkit import Chem

from ..basic import BuildingBlock
from ..basic import RetrosyntheticNode
from ..basic import RetrosyntheticPath
from ..basic.node import compute_inchi_key
from ..configs import BuildingBlockMatchPolicy
from .errors import UnsupportedBuildingBlockMatchPolicyError


class BuildingBlockChecker:
    """Checks whether all leaf molecules of a path are in a catalog.

    The catalog is loaded once from a plain-text / CSV SMILES file and
    stored as a frozen set for O(1) lookup. Keys are either RDKit-canonical
    SMILES (:attr:`~ursa.BuildingBlockMatchPolicy.SMILES`, default) or
    standard InChIKeys (:attr:`~ursa.BuildingBlockMatchPolicy.INCHI_KEY`)
    for tautomer-aware matching.

    :param catalog: Set of identity keys (canonical SMILES or InChIKeys).
    :type catalog: frozenset[str]
    :param match_policy: Which molecular identity to use for load/lookup.
    :type match_policy: BuildingBlockMatchPolicy
    """

    def __init__(
        self,
        catalog: frozenset[str],
        match_policy: BuildingBlockMatchPolicy | str = BuildingBlockMatchPolicy.SMILES,
    ) -> None:
        """Initialize BuildingBlockChecker with a pre-loaded catalog.

        Prefer :meth:`from_file` for constructing from a SMILES file.

        :param catalog: Frozen set of building-block identity keys.
        :type catalog: frozenset[str]
        :param match_policy: Catalog keying policy.
        :type match_policy: BuildingBlockMatchPolicy | str
        """
        self._catalog = catalog
        self._match_policy = BuildingBlockMatchPolicy(match_policy)

    @property
    def match_policy(self) -> BuildingBlockMatchPolicy:
        """Active building-block identity policy."""
        return self._match_policy

    @classmethod
    def from_file(
        cls,
        catalog_path: str | PathLike,
        match_policy: BuildingBlockMatchPolicy | str = BuildingBlockMatchPolicy.SMILES,
    ) -> BuildingBlockChecker:
        """Construct a :class:`BuildingBlockChecker` from a SMILES catalog.

        Each non-empty, non-comment line (or CSV ``smiles`` cell) is parsed
        and stored under the key dictated by ``match_policy``.

        :param catalog_path: Path to the building-block SMILES catalog file.
        :type catalog_path: str | PathLike
        :param match_policy: Catalog keying policy (default: canonical SMILES).
        :type match_policy: BuildingBlockMatchPolicy | str

        :return: A new :class:`BuildingBlockChecker` instance.
        :rtype: BuildingBlockChecker
        """
        policy = BuildingBlockMatchPolicy(match_policy)
        obj = cls(frozenset(), match_policy=policy)
        obj._catalog = obj._load_catalog(catalog_path)
        return obj

    def check(self, path: RetrosyntheticPath) -> tuple[BuildingBlock, ...]:
        """Check all starting-material (leaf) nodes of ``path`` against the catalog.

        Collects all nodes where
        :attr:`~ursa.RetrosyntheticNode.is_starting_material` is ``True``
        via :meth:`~ursa.RetrosyntheticPath.get_starting_materials` and looks
        up each molecule under the active :attr:`match_policy`.

        :param path: The retrosynthetic path whose starting materials to check.
        :type path: RetrosyntheticPath

        :return: Tuple of :class:`~ursa.BuildingBlock` objects, one per
            starting-material node, with ``found_in_catalog`` set accordingly.
        :rtype: tuple[BuildingBlock, ...]
        """
        return tuple(
            BuildingBlock(
                smiles=node.smiles,
                found_in_catalog=self._key_in_catalog(self._node_key(node)),
            )
            for node in path.get_starting_materials()
        )

    def _key_in_catalog(self, key: str) -> bool:
        """Return ``True`` when a non-empty ``key`` is in the catalog."""
        return bool(key) and key in self._catalog

    def _node_key(self, node: RetrosyntheticNode) -> str:
        """Return the catalog key for a leaf node.

        :raises UnsupportedBuildingBlockMatchPolicyError: If
            ``self._match_policy`` is not a known policy.
        """
        if self._match_policy is BuildingBlockMatchPolicy.INCHI_KEY:
            return node.inchi_key
        if self._match_policy is BuildingBlockMatchPolicy.SMILES:
            return node.canonical_smiles
        raise UnsupportedBuildingBlockMatchPolicyError(self._match_policy)

    def _smiles_key(self, smiles: str) -> str:
        """Return the catalog key for a raw catalog SMILES string.

        Uses uncached InChIKey conversion so a large catalog load does not
        pin every entry in the process-wide leaf cache.

        :raises UnsupportedBuildingBlockMatchPolicyError: If
            ``self._match_policy`` is not a known policy.
        """
        if self._match_policy is BuildingBlockMatchPolicy.INCHI_KEY:
            return compute_inchi_key(smiles)
        if self._match_policy is BuildingBlockMatchPolicy.SMILES:
            mol = Chem.MolFromSmiles(smiles)
            return Chem.MolToSmiles(mol) if mol is not None else ""
        raise UnsupportedBuildingBlockMatchPolicyError(self._match_policy)

    def _load_catalog(self, catalog_path: str | PathLike) -> frozenset[str]:
        """Load and return a frozen set of identity keys from a catalog file.

        Accepts two formats:
        - **CSV** with a ``smiles`` column header (e.g. ``bb_id,smiles``).
        - **Plain text** — one SMILES per line; skips blanks and ``#`` comments.

        Rows that fail SMILES parse (or InChIKey generation when that policy
        is active) are skipped.

        :param catalog_path: Path to the building-block catalog file.
        :type catalog_path: str | PathLike

        :return: Frozen set of identity keys.
        :rtype: frozenset[str]
        """
        keys: set[str] = set()
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
            key = self._smiles_key(smi)
            if key:
                keys.add(key)
        return frozenset(keys)
