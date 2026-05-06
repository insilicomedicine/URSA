import json
import logging

import pytest

from ursa.logging.formatter import JsonFormatter


@pytest.fixture
def formatter() -> JsonFormatter:
    return JsonFormatter()


def make_record(
    name: str = "test.logger",
    level: int = logging.INFO,
    msg: str = "hello",
    payload: dict | None = None,
) -> logging.LogRecord:
    record = logging.LogRecord(
        name=name,
        level=level,
        pathname="",
        lineno=0,
        msg=msg,
        args=(),
        exc_info=None,
    )
    if payload is not None:
        record._payload = payload  # type: ignore[attr-defined]
    return record


class TestFormat:
    def test_output_is_valid_json(self, formatter):
        record = make_record()
        result = formatter.format(record)
        assert json.loads(result)

    def test_timestamp_is_present(self, formatter):
        record = make_record()
        data = json.loads(formatter.format(record))
        assert "timestamp" in data

    def test_timestamp_is_iso8601(self, formatter):
        from datetime import datetime

        record = make_record()
        data = json.loads(formatter.format(record))
        # Should parse without error
        datetime.fromisoformat(data["timestamp"])

    def test_timestamp_matches_record_created(self, formatter):
        from datetime import datetime, timezone

        record = make_record()
        data = json.loads(formatter.format(record))
        expected = datetime.fromtimestamp(record.created, tz=timezone.utc).isoformat()
        assert data["timestamp"] == expected

    def test_level_field(self, formatter):
        record = make_record(level=logging.WARNING)
        data = json.loads(formatter.format(record))
        assert data["level"] == "WARNING"

    def test_logger_field(self, formatter):
        record = make_record(name="ursa.scoring")
        data = json.loads(formatter.format(record))
        assert data["logger"] == "ursa.scoring"

    def test_message_field(self, formatter):
        record = make_record(msg="test message")
        data = json.loads(formatter.format(record))
        assert data["message"] == "test message"

    def test_payload_fields_included(self, formatter):
        record = make_record(payload={"path_id": "p1", "score": 3.0})
        data = json.loads(formatter.format(record))
        assert data["path_id"] == "p1"
        assert data["score"] == pytest.approx(3.0)

    def test_missing_payload_produces_valid_json(self, formatter):
        record = make_record()  # no payload
        data = json.loads(formatter.format(record))
        assert "timestamp" in data
        assert "level" in data

    def test_multiple_records_are_independent(self, formatter):
        r1 = make_record(payload={"x": 1})
        r2 = make_record(payload={"y": 2})
        d1 = json.loads(formatter.format(r1))
        d2 = json.loads(formatter.format(r2))
        assert "x" in d1 and "y" not in d1
        assert "y" in d2 and "x" not in d2
