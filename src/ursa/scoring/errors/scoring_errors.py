from ...errors import UrsaError


class ScoringError(UrsaError):
    """Base class for reaction-scoring errors."""


class ReactionScoringError(ScoringError):
    """Raised when a reaction scorer fails to score a single reaction.

    :param reaction_smiles: The reaction SMILES that could not be scored.
    :type reaction_smiles: str
    :param msg: Human-readable error message.
    :type msg: str | None
    """

    def __init__(self, reaction_smiles: str, msg: str | None = None) -> None:
        """Initialize ReactionScoringError.

        :param reaction_smiles: The reaction SMILES that failed scoring.
        :type reaction_smiles: str
        :param msg: Human-readable error message.
        :type msg: str | None
        """
        self.reaction_smiles = reaction_smiles
        self.msg = msg or f"Failed to score reaction: {reaction_smiles!r}"
        super().__init__(self.msg)

    def __str__(self) -> str:
        return self.msg

    def __repr__(self) -> str:
        return (
            f"ReactionScoringError(reaction_smiles={self.reaction_smiles!r}, "
            f"msg={self.msg!r})"
        )


class EmptyVariantsError(ScoringError):
    """Raised when :class:`~ursa.BestPathSelector` receives no variants.

    At least one variant (the original path) must be present for
    selection to proceed.
    """

    def __init__(self, path_id: str) -> None:
        """Initialize EmptyVariantsError.

        :param path_id: Identifier of the path with no variants.
        :type path_id: str
        """
        self.path_id = path_id
        super().__init__(f"No variants available for path {path_id!r}")

    def __repr__(self) -> str:
        return f"EmptyVariantsError(path_id={self.path_id!r})"
