import json
import logging
from datetime import datetime
from datetime import timezone


class JsonFormatter(logging.Formatter):
    """JSON-lines formatter for structured log output.

    Converts each :class:`logging.LogRecord` into a single JSON object
    written on one line. The output is machine-readable and suitable for
    downstream parsing with tools such as ``jq`` or ``pandas``.

    Standard fields included in every record:

    * ``timestamp`` — ISO-8601 UTC time of the log event.
    * ``level`` — log level name (e.g. ``"INFO"``, ``"ERROR"``).
    * ``logger`` — name of the logger that emitted the record.
    * ``message`` — the formatted log message.

    Additional domain fields (e.g. ``path_id``, ``score``) are injected
    by :class:`~ursa.logging.Logger` via the record's ``extra`` dict and
    are included verbatim in the output JSON object.
    """

    def format(self, record: logging.LogRecord) -> str:
        """Serialize ``record`` to a JSON-lines string.

        :param record: The log record to format.
        :type record: logging.LogRecord

        :return: A single-line JSON string terminated by ``\\n``.
        :rtype: str
        """
        payload: dict = {
            "timestamp": datetime.fromtimestamp(
                record.created, tz=timezone.utc
            ).isoformat(),
            "level": record.levelname,
            "logger": record.name,
            "message": record.getMessage(),
        }
        payload.update(getattr(record, "_payload", {}))
        return json.dumps(payload)
