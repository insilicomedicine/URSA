from enum import Enum


class PathScoringConfig(Enum):
    """Configuration for scoring retrosynthetic paths.

    ``pass_threshold`` defines the minimum reaction-scorer value for a
    step to be considered as passed. Steps with a score strictly below
    this value are counted as failed when selecting the best collapsed
    variant and when computing ``is_route_solved``.

    ``failed_step_score`` is the sentinel value returned by the scorer
    for reactions that could not be processed (e.g. mapping failures).
    """

    pass_threshold = 1.0
    failed_step_score = 0.0
