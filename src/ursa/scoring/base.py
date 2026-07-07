from typing import Protocol


class ReactionScore(Protocol):
    """Dual-variant score for a single reaction.

    Mirrors ``chemcensor.ScoreResult``: the same reaction is scored both
    with and without functional-group matching in one evaluation. URSA
    uses the functional-group-agnostic score for Solv-1 and the
    functional-group-aware score for Solv-2.

    :param with_functional_groups: Score requiring each reaction center's
        functional-group sub-signature to match the reference (Solv-2).
    :type with_functional_groups: float
    :param without_functional_groups: Score based on reaction-center
        presence only, ignoring functional groups (Solv-1).
    :type without_functional_groups: float
    """

    with_functional_groups: float
    without_functional_groups: float


class ReactionScorer(Protocol):
    """Protocol for single-reaction scorers.

    Any object implementing this protocol can be used as a reaction
    scorer inside :class:`~ursa.PathScorer`. The canonical implementation
    is :class:`chemcensor.ChemCensor`.
    """

    def evaluate(self, reaction_smiles: str) -> ReactionScore:
        """Evaluate a reaction, returning both functional-group variants.

        :param reaction_smiles: Forward reaction SMILES to score.
        :type reaction_smiles: str

        :return: Object exposing ``with_functional_groups`` and
            ``without_functional_groups`` float scores. Higher values
            indicate better-known reactions; negative values (e.g.
            ``-1.0``) indicate scoring failures.
        :rtype: ReactionScore
        """
        ...
