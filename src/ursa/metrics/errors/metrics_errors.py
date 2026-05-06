from ...errors import UrsaError


class MetricsError(UrsaError):
    """Base class for dataset-metrics errors."""


class EmptyDatasetError(MetricsError):
    """Raised when metrics are requested for an empty dataset.

    :class:`~ursa.DatasetMetricsCalculator` requires at least one
    :class:`~ursa.PathResult` to compute meaningful metrics.
    """

    def __init__(self) -> None:
        """Initialize EmptyDatasetError."""
        super().__init__("Cannot compute metrics for an empty dataset")

    def __repr__(self) -> str:
        return "EmptyDatasetError()"
