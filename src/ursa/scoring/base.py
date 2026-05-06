from typing import Protocol


class ReactionScorer(Protocol):
    """Protocol for single-reaction scorers.

    Any object implementing this protocol can be used as a reaction
    scorer inside :class:`~ursa.PathScorer`. The canonical implementation
    is :class:`chemcensor.ChemCensor`, but any callable object with a
    matching ``score`` signature is accepted.
    """

    def score(self, reaction_smiles: str) -> float:
        """Score a single reaction.

        :param reaction_smiles: Forward reaction SMILES to score.
        :type reaction_smiles: str

        :return: Numeric score. Higher values indicate better-known
            reactions. Negative values (e.g. ``-1.0``) indicate scoring
            failures.
        :rtype: float
        """
        ...
