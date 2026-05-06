from enum import Enum


class ValidationConfig(Enum):
    """Configuration for path validation and collapsing.

    ``min_path_depth`` is the minimum acceptable depth of a
    retrosynthetic tree. Paths with fewer steps than this value are
    rejected by :class:`~ursa.PathConsistencyChecker`.

    ``combination_step_limit`` is the maximum number of steps for which
    :class:`~ursa.PathCollapser` will enumerate all valid collapse
    combinations. Paths exceeding this limit are returned as-is
    (only the original tree, without collapsed variants) to avoid
    combinatorial explosion.
    """

    min_path_depth = 1
    combination_step_limit = 20
