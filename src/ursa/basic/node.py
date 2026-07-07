from __future__ import annotations

from dataclasses import dataclass
from dataclasses import field
from functools import cached_property
from functools import lru_cache

from rdkit import Chem


@lru_cache(maxsize=None)
def _canonicalize(smiles: str) -> tuple[str, bool]:
    """Return ``(canonical_smiles, is_valid)`` for ``smiles``, memoised.

    Canonicalisation is a pure function of the input string, so results
    are cached: :class:`~ursa.PathCollapser` rebuilds the same molecules
    many times across collapsed variants, and re-parsing each one with
    RDKit dominates the prepare phase otherwise.

    :param smiles: Input SMILES string.
    :type smiles: str
    :return: RDKit-canonical SMILES (``""`` if unparsable) and a validity
        flag.
    :rtype: tuple[str, bool]
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return "", False
    return Chem.MolToSmiles(mol), True


@dataclass(frozen=True)
class RetrosyntheticNode:
    """A node in a retrosynthetic tree.

    Each node represents a molecule that is either disconnected into
    precursors (internal node) or is a terminal building block (leaf).

    :param smiles: SMILES string of the molecule at this node.
    :type smiles: str
    :param children: Precursor nodes produced by the retrosynthetic step.
        Empty tuple indicates a starting-material (leaf) node.
    :type children: tuple[RetrosyntheticNode, ...]
    :param canonical_smiles: RDKit-canonical SMILES; empty string if the
        molecule could not be parsed.
    :type canonical_smiles: str
    :param is_valid: ``True`` if RDKit could parse ``smiles``.
    :type is_valid: bool
    """

    smiles: str
    children: tuple[RetrosyntheticNode, ...] = field(default_factory=tuple)
    canonical_smiles: str = field(init=False)
    is_valid: bool = field(init=False)

    def __post_init__(self) -> None:
        canonical_smiles, is_valid = _canonicalize(self.smiles)
        object.__setattr__(self, "canonical_smiles", canonical_smiles)
        object.__setattr__(self, "is_valid", is_valid)

    @classmethod
    def _from_mol_dict(cls, mol_node: dict) -> RetrosyntheticNode:
        """Construct a node recursively from a RetroCast molecule dict.

        Expected ``mol_node`` structure (RetroCast schema 2)::

            {
                "smiles": "<SMILES>",
                "inchikey": "<InChIKey>",
                "product_of": {
                    "reactants": [ { ... same structure ... }, ... ]
                } | null
            }

        A molecule whose ``product_of`` is ``null`` or has empty
        ``reactants`` becomes a leaf node with no children.

        :param mol_node: Molecule object from a RetroCast JSON route.
        :type mol_node: dict

        :return: Fully populated :class:`RetrosyntheticNode` for this molecule
            and all its descendants.
        :rtype: RetrosyntheticNode
        """
        smiles = mol_node["smiles"]
        step = mol_node.get("product_of")
        if not step:
            return cls(smiles=smiles)
        reactants = step.get("reactants") or ()
        if not reactants:
            return cls(smiles=smiles)
        children = tuple(cls._from_mol_dict(r) for r in reactants)
        return cls(smiles=smiles, children=children)

    @cached_property
    def reaction_smiles(self) -> str:
        """Forward reaction SMILES derived from children and this node's molecule.

        Format: ``"reactant1.reactant2>>product"``. Returns an empty string
        for leaf nodes (starting materials with no children).

        :return: Unmapped forward reaction SMILES, or ``""`` for leaf nodes.
        :rtype: str
        """
        if not self.children:
            return ""
        reactants = ".".join(c.canonical_smiles for c in self.children)
        return f"{reactants}>>{self.canonical_smiles}"

    @cached_property
    def is_starting_material(self) -> bool:
        """Return ``True`` if this node is a leaf (starting material).

        A node is a starting material when it has no children, i.e. it is
        not further disconnected in the retrosynthetic tree. This matches
        the ``is_starting_material`` flag used in the Syntheseus output format.

        Whether a starting material is commercially available (i.e. a
        *building block*) is determined separately by
        :class:`~ursa.BuildingBlockChecker` and stored in
        :class:`~ursa.BuildingBlock`.

        :return: ``True`` if the node is a starting material, ``False`` otherwise.
        :rtype: bool
        """
        return len(self.children) == 0

    def __hash__(self) -> int:
        cached = self.__dict__.get("_hash")
        if cached is not None:
            return cached
        value = hash((self.canonical_smiles, self.children))
        object.__setattr__(self, "_hash", value)
        return value

    def __eq__(self, other: object) -> bool:
        if not isinstance(other, RetrosyntheticNode):
            return NotImplemented
        return (
            self.canonical_smiles == other.canonical_smiles
            and self.children == other.children
        )
