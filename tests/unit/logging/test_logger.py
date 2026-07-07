import logging
from unittest.mock import MagicMock
from unittest.mock import patch

import pytest

from tests.conftest import make_variant
from ursa.logging.formatter import JsonFormatter
from ursa.logging.logger import Logger


# ── helpers ──────────────────────────────────────────────────────────────────


def _payload(mock_std_logger: MagicMock) -> dict:
    """Return the _payload dict from the most recent log() call."""
    return mock_std_logger.log.call_args.kwargs["extra"]["_payload"]


@pytest.fixture
def mock_std_logger():
    return MagicMock()


@pytest.fixture
def logger(mock_std_logger):
    with patch("logging.getLogger", return_value=mock_std_logger):
        return Logger("test")


# ── LogStep ───────────────────────────────────────────────────────────────────


class TestLogStep:
    def test_calls_underlying_log(self, logger, mock_std_logger, path_1step):
        vr = make_variant(path_1step, (3.0,))
        logger.log_step("p1", vr.step_results[0])
        mock_std_logger.log.assert_called_once()

    def test_payload_path_id(self, logger, mock_std_logger, path_1step):
        vr = make_variant(path_1step, (3.0,))
        logger.log_step("p1", vr.step_results[0])
        assert _payload(mock_std_logger)["path_id"] == "p1"

    def test_payload_smiles(self, logger, mock_std_logger, path_1step):
        vr = make_variant(path_1step, (3.0,))
        sr = vr.step_results[0]
        logger.log_step("p1", sr)
        assert _payload(mock_std_logger)["smiles"] == sr.node.smiles

    def test_payload_reaction_smiles(self, logger, mock_std_logger, path_1step):
        vr = make_variant(path_1step, (3.0,))
        sr = vr.step_results[0]
        logger.log_step("p1", sr)
        assert _payload(mock_std_logger)["reaction_smiles"] == sr.node.reaction_smiles

    def test_payload_scores(self, logger, mock_std_logger, path_1step):
        vr = make_variant(path_1step, (3.0,))
        logger.log_step("p1", vr.step_results[0])
        p = _payload(mock_std_logger)
        assert p["score_without_fg"] == pytest.approx(3.0)
        assert p["score_with_fg"] == pytest.approx(3.0)

    def test_payload_passes(self, logger, mock_std_logger, path_1step):
        vr = make_variant(path_1step, (3.0,))
        logger.log_step("p1", vr.step_results[0])
        p = _payload(mock_std_logger)
        assert p["passes_solv_1"] is True
        assert p["passes_solv_2"] is True


# ── LogVariant ────────────────────────────────────────────────────────────────


class TestLogVariant:
    def test_payload_num_steps(self, logger, mock_std_logger, path_3step_linear):
        vr = make_variant(path_3step_linear, (3.0, 2.0, 1.0))
        logger.log_variant("p3", vr)
        assert _payload(mock_std_logger)["num_steps"] == 3

    def test_payload_mean_scores(self, logger, mock_std_logger, path_3step_linear):
        vr = make_variant(path_3step_linear, (3.0, 1.0, 2.0))
        logger.log_variant("p3", vr)
        p = _payload(mock_std_logger)
        assert p["mean_score_without_fg"] == pytest.approx(2.0)
        assert p["mean_score_with_fg"] == pytest.approx(2.0)

    def test_payload_all_steps_pass(self, logger, mock_std_logger, path_1step):
        vr = make_variant(path_1step, (3.0,))
        logger.log_variant("p1", vr)
        p = _payload(mock_std_logger)
        assert p["all_steps_pass_solv_1"] is True
        assert p["all_steps_pass_solv_2"] is True


# ── LogPath ───────────────────────────────────────────────────────────────────


class TestLogPath:
    def _make_path_result(self, path, scores):
        from tests.conftest import make_path_result

        return make_path_result(path, scores)

    def test_payload_path_id(self, logger, mock_std_logger, path_1step):
        pr = self._make_path_result(path_1step, (3.0,))
        logger.log_path(pr)
        assert _payload(mock_std_logger)["path_id"] == path_1step.path_id

    def test_payload_is_consistent(self, logger, mock_std_logger, path_1step):
        pr = self._make_path_result(path_1step, (3.0,))
        logger.log_path(pr)
        assert _payload(mock_std_logger)["is_consistent"] is True

    def test_payload_all_bb_found(self, logger, mock_std_logger, path_1step):
        pr = self._make_path_result(path_1step, (3.0,))
        logger.log_path(pr)
        assert _payload(mock_std_logger)["all_bb_found"] is True

    def test_payload_passes_solv_levels(self, logger, mock_std_logger, path_1step):
        pr = self._make_path_result(path_1step, (3.0,))
        logger.log_path(pr)
        p = _payload(mock_std_logger)
        assert p["passes_solv_0"] is True
        assert p["passes_solv_1"] is True
        assert p["passes_solv_2"] is True

    def test_payload_best_variant_summary(self, logger, mock_std_logger, path_1step):
        pr = self._make_path_result(path_1step, (3.0,))
        logger.log_path(pr)
        p = _payload(mock_std_logger)
        assert "num_steps" in p
        assert "is_no_synthesis" in p
        assert "mean_score_with_fg" in p


# ── LogDataset ────────────────────────────────────────────────────────────────


class TestLogDataset:
    def _make_dataset_result(self, path, scores):
        from tests.conftest import make_path_result
        from ursa.basic.result import DatasetMetrics, DatasetResult

        pr = make_path_result(path, scores)
        mean = sum(scores) / len(scores)
        metrics = DatasetMetrics(
            total_molecules=1,
            molecules_with_route=1,
            routes_no_synthesis=0,
            routes_solv_0=1,
            routes_solv_1=1,
            routes_solv_2=1,
            solv_0=1.0,
            solv_1=1.0,
            solv_2=1.0,
            mean_score_without_fg=mean,
            mean_score_with_fg=mean,
        )
        return DatasetResult(path_results=(pr,), metrics=metrics)

    def test_payload_contains_all_metric_fields(
        self, logger, mock_std_logger, path_1step
    ):
        dr = self._make_dataset_result(path_1step, (3.0,))
        logger.log_dataset(dr)
        p = _payload(mock_std_logger)
        expected_keys = {
            "total_molecules",
            "molecules_with_route",
            "routes_no_synthesis",
            "routes_solv_0",
            "routes_solv_1",
            "routes_solv_2",
            "solv_0",
            "solv_1",
            "solv_2",
            "mean_score_without_fg",
            "mean_score_with_fg",
        }
        assert expected_keys <= p.keys()

    def test_payload_values_match_metrics(self, logger, mock_std_logger, path_1step):
        dr = self._make_dataset_result(path_1step, (3.0,))
        logger.log_dataset(dr)
        p = _payload(mock_std_logger)
        m = dr.metrics
        assert p["total_molecules"] == m.total_molecules
        assert p["routes_solv_2"] == m.routes_solv_2
        assert p["solv_2"] == pytest.approx(m.solv_2)


# ── GetLogger ─────────────────────────────────────────────────────────────────


class TestGetLogger:
    def test_returns_logger_instance(self):
        from ursa.logging import get_logger

        result = get_logger("test.get_logger")
        assert isinstance(result, Logger)

    def test_handler_is_attached(self):
        from ursa.logging import get_logger

        result = get_logger("test.handler_attached")
        assert len(result._logger.handlers) >= 1

    def test_json_formatter_used_by_default(self):
        from ursa.logging import get_logger

        result = get_logger("test.json_fmt")
        handlers = result._logger.handlers
        assert any(isinstance(h.formatter, JsonFormatter) for h in handlers)

    def test_level_is_set(self):
        from ursa.logging import get_logger

        result = get_logger("test.level_set")
        assert result._logger.level == logging.INFO
