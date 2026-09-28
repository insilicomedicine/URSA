from __future__ import annotations

from rdkit import Chem


class Molecule:
    """Minimal molecule wrapper exposing a canonical RDKit mol for drawing.

    The CDXML renderer only needs a canonical, atom-map-free RDKit molecule to
    lay out atoms and bonds. This is a lightweight stand-in for the richer
    ``Molecule`` class used by the legacy benchmark, providing just the
    :attr:`canonical_rdmol` property the renderer relies on.

    :param smiles: SMILES string of the molecule.
    :type smiles: str
    """

    def __init__(self, smiles: str) -> None:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            raise ValueError(f"Invalid SMILES: {smiles!r}")
        canonical = Chem.MolToSmiles(mol, ignoreAtomMapNumbers=True)
        self._canonical_mol = Chem.MolFromSmiles(canonical)
        if self._canonical_mol is None:
            raise ValueError(f"Could not canonicalize SMILES: {smiles!r}")

    @property
    def canonical_rdmol(self) -> Chem.Mol:
        """Canonical RDKit molecule without atom maps.

        :return: The canonical RDKit molecule.
        :rtype: rdkit.Chem.Mol
        """
        return self._canonical_mol
