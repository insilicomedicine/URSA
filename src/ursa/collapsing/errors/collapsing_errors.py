from ...errors import UrsaError


class CollapsingError(UrsaError):
    """Base class for path-collapsing errors."""


class CyclicTreeError(CollapsingError):
    """Raised when the retrosynthetic tree contains a cycle.

    A valid retrosynthetic tree is acyclic. A cycle would cause
    :class:`~ursa.PathCollapser` to loop indefinitely during traversal.
    """

    def __init__(self, path_id: str) -> None:
        """Initialize CyclicTreeError.

        :param path_id: Identifier of the path with the cycle.
        :type path_id: str
        """
        self.path_id = path_id
        super().__init__(f"Path {path_id!r} contains a cycle and cannot be collapsed")

    def __repr__(self) -> str:
        return f"CyclicTreeError(path_id={self.path_id!r})"
