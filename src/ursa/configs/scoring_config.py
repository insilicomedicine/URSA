from enum import Enum


class BestPathSelectionPolicy(str, Enum):
    """How :class:`~ursa.BestPathSelector` ranks variants after failed steps.

    Failed-step count is always the primary key (Solv-N existence check).
    The policy only reorders the remaining tie-breakers:

    * ``MEAN_SCORE`` (default) — highest mean ChemCensor score, then
      shortest path.
    * ``PATH_LENGTH`` — shortest path, then highest mean ChemCensor score.
    """

    MEAN_SCORE = "mean_score"
    PATH_LENGTH = "path_length"


class PathScoringConfig(Enum):
    """Configuration for scoring retrosynthetic paths.

    ``failed_step_score`` is the normalised floor for a reaction score:
    negative scorer values (processing failures) are clamped to it.

    ``pass_threshold`` is the minimum (normalised) reaction score a step
    must **strictly exceed** to pass a Solv level. Steps scoring at or
    below it are counted as failed both when deciding ``passes_solv_1/2``
    and when the best-variant selector ranks variants. The default
    (``0.0``) matches the Solv-N definition "every step scores > 0"; raise
    it to demand a higher per-step confidence.

    Variant ranking after the failed-step key is controlled by
    :class:`BestPathSelectionPolicy` (passed to
    :class:`~ursa.BestPathSelector` / :class:`~ursa.Ursa`), not by this
    enum.
    """

    failed_step_score = 0.0
    pass_threshold = 0.0
