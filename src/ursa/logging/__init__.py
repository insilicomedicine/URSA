import logging as _logging

from ..configs import LoggingConfig
from .formatter import JsonFormatter
from .logger import Logger


def get_logger(
    name: str, config: LoggingConfig = LoggingConfig.default_level
) -> Logger:
    """Return a configured :class:`Logger` instance.

    Attaches a :class:`JsonFormatter` handler to the underlying
    :class:`logging.Logger` when :attr:`LoggingConfig.json_format` is
    ``True``, otherwise falls back to a plain-text formatter.

    :param name: Logger name (typically ``__name__`` of the caller).
    :type name: str
    :param config: Logging configuration controlling level and format.
        Defaults to :attr:`LoggingConfig.default_level`.
    :type config: LoggingConfig

    :return: Configured :class:`Logger` instance.
    :rtype: Logger
    """
    logger = Logger(name, config)
    handler = _logging.StreamHandler()
    if LoggingConfig.json_format.value:
        handler.setFormatter(JsonFormatter())
    else:
        handler.setFormatter(_logging.Formatter())
    logger._logger.addHandler(handler)
    logger._logger.setLevel(config.value)
    return logger


__all__ = [
    "Logger",
    "JsonFormatter",
    "get_logger",
]
