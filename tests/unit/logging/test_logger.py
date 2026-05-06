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

    def test_payload_score(self, logger, mock_std_logger, path_1step):
        vr = make_variant(path_1step, (3.0,))
        logger.log_step("p1", vr.step_results[0])
        assert _payload(mock_std_logger)["score"] == pytest.approx(3.0)

    def test_payload_passed(self, logger, mock_std_logger, path_1step):
        vr = make_variant(path_1step, (3.0,))
        logger.log_step("p1", vr.step_results[0])
        assert _payload(mock_std_logger)["passed"] is True


# ── LogVariant ────────────────────────────────────────────────────────────────


class TestLogVariant:
    def test_payload_num_steps(self, logger, mock_std_logger, path_3step_linear):
        vr = make_variant(path_3step_linear, (3.0, 2.0, 1.0))
        logger.log_variant("p3", vr)
        assert _payload(mock_std_logger)["num_steps"] == 3

    def test_payload_percent_found(self, logger, mock_std_logger, path_1step):
        vr = make_variant(path_1step, (3.0,))
        logger.log_variant("p1", vr)
        assert _payload(mock_std_logger)["percent_found"] == pytest.approx(1.0)

    def test_payload_chemcensor_per_route(
        self, logger, mock_std_logger, path_3step_linear
    ):
        vr = make_variant(path_3step_linear, (3.0, 1.0, 2.0))
        logger.log_variant("p3", vr)
        assert _payload(mock_std_logger)["chemcensor_per_route"] == pytest.approx(2.0)

    def test_payload_all_steps_passed(self, logger, mock_std_logger, path_1step):
        vr = make_variant(path_1step, (3.0,))
        logger.log_variant("p1", vr)
        assert _payload(mock_std_logger)["all_steps_passed"] is True


# ── LogPath ───────────────────────────────────────────────────────────────────


class TestLogPath:
    def _make_path_result(self, path, scores):
        from ursa.basic.building_block import BuildingBlock
        from ursa.basic.result import PathResult

        vr = make_variant(path, scores)
        return PathResult(
            original_path=path,
            is_consistent=True,
            starting_materials=(BuildingBlock(smiles="CC", found_in_catalog=True),),
            all_bb_found=True,
            best_variant=vr,
            is_route_solved=True,
        )

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

    def test_payload_is_route_solved(self, logger, mock_std_logger, path_1step):
        pr = self._make_path_result(path_1step, (3.0,))
        logger.log_path(pr)
        assert _payload(mock_std_logger)["is_route_solved"] is True

    def test_payload_best_variant_summary(self, logger, mock_std_logger, path_1step):
        pr = self._make_path_result(path_1step, (3.0,))
        logger.log_path(pr)
        p = _payload(mock_std_logger)
        assert "num_steps" in p
        assert "percent_found" in p
        assert "chemcensor_per_route" in p
        assert "all_steps_passed" in p


# ── LogDataset ────────────────────────────────────────────────────────────────


class TestLogDataset:
    def _make_dataset_result(self, path, scores):
        from ursa.basic.building_block import BuildingBlock
        from ursa.basic.result import DatasetMetrics, DatasetResult, PathResult

        vr = make_variant(path, scores)
        pr = PathResult(
            original_path=path,
            is_consistent=True,
            starting_materials=(BuildingBlock(smiles="CC", found_in_catalog=True),),
            all_bb_found=True,
            best_variant=vr,
            is_route_solved=True,
        )
        metrics = DatasetMetrics(
            total_molecules=1,
            molecules_with_route=1,
            solved_routes=1,
            passed_steps=len(scores),
            total_steps=len(scores),
            solv_2=1.0,
            mean_chemcensor_score=sum(scores) / len(scores),
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
            "solved_routes",
            "passed_steps",
            "total_steps",
            "solv_2",
            "mean_chemcensor_score",
        }
        assert expected_keys <= p.keys()

    def test_payload_values_match_metrics(self, logger, mock_std_logger, path_1step):
        dr = self._make_dataset_result(path_1step, (3.0,))
        logger.log_dataset(dr)
        p = _payload(mock_std_logger)
        m = dr.metrics
        assert p["total_molecules"] == m.total_molecules
        assert p["solved_routes"] == m.solved_routes
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
