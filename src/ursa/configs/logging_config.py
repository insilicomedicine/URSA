from enum import Enum


class LoggingConfig(Enum):
    """Configuration for the structured logger.

    ``default_level`` is the logging level used when no explicit level
    is passed to :func:`~ursa.logging.get_logger`.

    ``json_format`` controls whether log records are emitted as
    machine-readable JSON lines (``True``) or as plain human-readable
    text (``False``).
    """

    default_level = "INFO"
    json_format = True
