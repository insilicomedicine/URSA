from ...errors import UrsaError


class PathError(UrsaError):
    """Base class for retrosynthetic path errors."""


class EmptyPathError(PathError):
    """Raised when a retrosynthetic path contains no steps."""

    def __init__(self, path_id: str) -> None:
        """Initialize EmptyPathError.

        :param path_id: Identifier of the empty path.
        :type path_id: str
        """
        self.path_id = path_id
        super().__init__(f"Path {path_id!r} contains no steps")

    def __repr__(self) -> str:
        return f"EmptyPathError(path_id={self.path_id!r})"


class InvalidNodeError(PathError):
    """Raised when a retrosynthetic node has an invalid structure."""

    def __init__(self, smiles: str, msg: str | None = None) -> None:
        """Initialize InvalidNodeError.

        :param smiles: SMILES of the molecule at the invalid node.
        :type smiles: str
        :param msg: Human-readable error message.
        :type msg: str | None
        """
        self.smiles = smiles
        self.msg = msg or f"Invalid node for molecule {smiles!r}"
        super().__init__(self.msg)

    def __str__(self) -> str:
        return self.msg

    def __repr__(self) -> str:
        return f"InvalidNodeError(smiles={self.smiles!r}, " f"msg={self.msg!r})"
