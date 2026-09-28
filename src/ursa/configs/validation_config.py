from enum import Enum


class BuildingBlockMatchPolicy(str, Enum):
    """How :class:`~ursa.BuildingBlockChecker` keys the BB catalog.

    * ``SMILES`` (default) — RDKit-canonical SMILES (legacy behaviour).
    * ``INCHI_KEY`` — standard InChIKey; tautomer SMILES of the same
      building block match.
    """

    SMILES = "smiles"
    INCHI_KEY = "inchi_key"


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

    Catalog identity for stock termination is controlled by
    :class:`BuildingBlockMatchPolicy` (passed to
    :class:`~ursa.BuildingBlockChecker` / :class:`~ursa.Ursa`).
    """

    min_path_depth = 1
    combination_step_limit = 20
