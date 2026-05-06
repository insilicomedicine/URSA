from .validation_errors import ValidationError


class ConsistencyError(ValidationError):
    """Base class for path-consistency errors."""


class DisconnectedTreeError(ConsistencyError):
    """Raised when the retrosynthetic tree is not fully connected.

    This happens when some molecules in the tree are unreachable from
    the root via the product-to-reactant relationships.
    """

    def __init__(self, path_id: str, unreachable_smiles: list[str]) -> None:
        """Initialize DisconnectedTreeError.

        :param path_id: Identifier of the inconsistent path.
        :type path_id: str
        :param unreachable_smiles: SMILES of molecules that could not
            be reached from the root.
        :type unreachable_smiles: list[str]
        """
        self.path_id = path_id
        self.unreachable_smiles = unreachable_smiles
        super().__init__(
            f"Path {path_id!r} has {len(unreachable_smiles)} unreachable "
            f"molecule(s): {unreachable_smiles[:3]}"
        )

    def __repr__(self) -> str:
        return (
            f"DisconnectedTreeError(path_id={self.path_id!r}, "
            f"unreachable_smiles={self.unreachable_smiles!r})"
        )


class MultipleRootsError(ConsistencyError):
    """Raised when the retrosynthetic tree has more than one root.

    A valid tree has exactly one molecule that is a product but never
    a reactant (the target molecule).
    """

    def __init__(self, path_id: str, root_smiles: list[str]) -> None:
        """Initialize MultipleRootsError.

        :param path_id: Identifier of the inconsistent path.
        :type path_id: str
        :param root_smiles: SMILES of all candidate root molecules.
        :type root_smiles: list[str]
        """
        self.path_id = path_id
        self.root_smiles = root_smiles
        super().__init__(
            f"Path {path_id!r} has {len(root_smiles)} root(s): {root_smiles}"
        )

    def __repr__(self) -> str:
        return (
            f"MultipleRootsError(path_id={self.path_id!r}, "
            f"root_smiles={self.root_smiles!r})"
        )
