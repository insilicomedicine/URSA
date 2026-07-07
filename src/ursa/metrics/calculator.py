from ..basic import DatasetMetrics
from ..basic import PathResult


class DatasetMetricsCalculator:
    """Computes aggregate Solv-0/1/2 metrics over a dataset of evaluated paths.

    Aggregates per-path :class:`~ursa.PathResult` objects into a single
    :class:`~ursa.DatasetMetrics`. The ``total_molecules`` parameter is
    the denominator for every ``solv_N`` rate, so targets for which no
    route was provided correctly lower the rates.
    """

    def calculate(
        self,
        path_results: tuple[PathResult, ...],
        total_molecules: int,
    ) -> DatasetMetrics:
        """Compute dataset-level metrics from ``path_results``.

        :param path_results: Per-path evaluation results to aggregate.
        :type path_results: tuple[PathResult, ...]
        :param total_molecules: Total number of target molecules in the
            input dataset, including those without a route. Denominator
            for every ``solv_N`` rate.
        :type total_molecules: int

        :return: Aggregated dataset metrics.
        :rtype: DatasetMetrics
        """
        molecules_with_route = len(path_results)
        routes_no_synthesis = sum(1 for pr in path_results if pr.is_no_synthesis)
        routes_solv_0 = sum(1 for pr in path_results if pr.passes_solv_0)
        routes_solv_1 = sum(1 for pr in path_results if pr.passes_solv_1)
        routes_solv_2 = sum(1 for pr in path_results if pr.passes_solv_2)

        def _rate(count: int) -> float:
            return count / total_molecules if total_molecules > 0 else 0.0

        steps_without_fg = [
            sr.score_without_fg
            for pr in path_results
            for sr in pr.best_variant_solv_1.step_results
        ]
        steps_with_fg = [
            sr.score_with_fg
            for pr in path_results
            for sr in pr.best_variant_solv_2.step_results
        ]
        mean_score_without_fg = (
            sum(steps_without_fg) / len(steps_without_fg) if steps_without_fg else 0.0
        )
        mean_score_with_fg = (
            sum(steps_with_fg) / len(steps_with_fg) if steps_with_fg else 0.0
        )

        return DatasetMetrics(
            total_molecules=total_molecules,
            molecules_with_route=molecules_with_route,
            routes_no_synthesis=routes_no_synthesis,
            routes_solv_0=routes_solv_0,
            routes_solv_1=routes_solv_1,
            routes_solv_2=routes_solv_2,
            solv_0=_rate(routes_solv_0),
            solv_1=_rate(routes_solv_1),
            solv_2=_rate(routes_solv_2),
            mean_score_without_fg=mean_score_without_fg,
            mean_score_with_fg=mean_score_with_fg,
        )
