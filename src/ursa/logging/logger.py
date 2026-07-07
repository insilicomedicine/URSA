import logging

from ..basic import DatasetResult
from ..basic import PathResult
from ..basic import StepResult
from ..basic import VariantResult
from ..configs import LoggingConfig


class Logger:
    """Structured logger for ursa evaluation events.

    Wraps Python's standard :mod:`logging` and emits one JSON-lines
    record per event via :class:`~ursa.logging.JsonFormatter`. Four
    granularity levels are supported:

    * **step** — one record per scored reaction.
    * **variant** — one record per collapsed path variant (after scoring).
    * **path** — one record per evaluated :class:`~ursa.RetrosyntheticPath`.
    * **dataset** — one summary record per :meth:`~ursa.Ursa.score_dataset`
      call.

    :param name: Logger name passed to :func:`logging.getLogger`.
    :type name: str
    :param config: Logging configuration controlling level and format.
    :type config: LoggingConfig
    """

    def __init__(
        self,
        name: str,
        config: LoggingConfig = LoggingConfig.default_level,
    ) -> None:
        """Initialize Logger.

        :param name: Logger name (typically ``__name__`` of the caller).
        :type name: str
        :param config: Logging configuration. Defaults to
            :attr:`LoggingConfig.default_level`.
        :type config: LoggingConfig
        """
        self._logger = logging.getLogger(name)
        self._config = config

    def log_step(self, path_id: str, result: StepResult) -> None:
        """Emit a structured log record for a single scored step.

        Includes ``path_id``, ``smiles``, ``reaction_smiles``, both score
        variants, and per-level pass flags.

        :param path_id: Identifier of the path this step belongs to.
        :type path_id: str
        :param result: Scoring result for the step.
        :type result: StepResult
        """
        self._emit(
            path_id=path_id,
            smiles=result.node.smiles,
            reaction_smiles=result.node.reaction_smiles,
            score_without_fg=result.score_without_fg,
            score_with_fg=result.score_with_fg,
            passes_solv_1=result.passes_solv_1,
            passes_solv_2=result.passes_solv_2,
        )

    def log_variant(self, path_id: str, result: VariantResult) -> None:
        """Emit a structured log record for a scored path variant.

        Includes ``path_id``, ``num_steps``, per-level pass flags, and
        mean scores.

        :param path_id: Identifier of the original path.
        :type path_id: str
        :param result: Scoring result for the variant.
        :type result: VariantResult
        """
        self._emit(
            path_id=path_id,
            num_steps=result.path.num_steps,
            all_steps_pass_solv_1=result.all_steps_pass_solv_1,
            all_steps_pass_solv_2=result.all_steps_pass_solv_2,
            mean_score_without_fg=result.mean_score_without_fg,
            mean_score_with_fg=result.mean_score_with_fg,
        )

    def log_path(self, result: PathResult) -> None:
        """Emit a structured log record for a fully evaluated path.

        Includes ``path_id``, ``is_consistent``, ``all_bb_found``,
        ``is_no_synthesis``, the three ``passes_solv_*`` flags, and summary
        scores from the Solv-2 best variant.

        :param result: Full evaluation result for the path.
        :type result: PathResult
        """
        bv = result.best_variant_solv_2
        self._emit(
            path_id=result.original_path.path_id,
            is_consistent=result.is_consistent,
            all_bb_found=result.all_bb_found,
            is_no_synthesis=result.is_no_synthesis,
            passes_solv_0=result.passes_solv_0,
            passes_solv_1=result.passes_solv_1,
            passes_solv_2=result.passes_solv_2,
            num_steps=bv.path.num_steps,
            mean_score_with_fg=bv.mean_score_with_fg,
        )

    def log_dataset(self, result: DatasetResult) -> None:
        """Emit a structured log record summarising dataset evaluation.

        Includes the Solv-0/1/2 rates and counts plus diagnostic mean
        scores from :class:`~ursa.DatasetMetrics`.

        :param result: Full dataset evaluation result.
        :type result: DatasetResult
        """
        m = result.metrics
        self._emit(
            total_molecules=m.total_molecules,
            molecules_with_route=m.molecules_with_route,
            routes_no_synthesis=m.routes_no_synthesis,
            routes_solv_0=m.routes_solv_0,
            routes_solv_1=m.routes_solv_1,
            routes_solv_2=m.routes_solv_2,
            solv_0=m.solv_0,
            solv_1=m.solv_1,
            solv_2=m.solv_2,
            mean_score_without_fg=m.mean_score_without_fg,
            mean_score_with_fg=m.mean_score_with_fg,
        )

    def _emit(self, **fields: object) -> None:
        """Forward ``fields`` to the underlying :mod:`logging` logger.

        Wraps the payload in the record's ``extra["_payload"]`` slot so
        that :class:`~ursa.logging.JsonFormatter` (or the CLI formatter)
        can serialise it without colliding with standard logging keys.

        :param fields: Arbitrary key-value payload to attach to the record.
        """
        level = logging.getLevelName(self._config.value)
        self._logger.log(level, "", extra={"_payload": fields})
