from ..basic import DatasetMetrics
from ..basic import PathResult


class DatasetMetricsCalculator:
    """Computes aggregate metrics over a dataset of evaluated paths.

    Aggregates per-path :class:`~ursa.PathResult` objects into a single
    :class:`~ursa.DatasetMetrics` instance. The ``total_molecules``
    parameter allows computing ``solv_2`` relative to the
    full input dataset (including molecules for which no route was
    provided).
    """

    def calculate(
        self,
        path_results: tuple[PathResult, ...],
        total_molecules: int,
    ) -> DatasetMetrics:
        """Compute dataset-level metrics from ``path_results``.

        Counts solved routes, passed/total steps and derives
        ``solv_2`` and ``mean_chemcensor_score``.

        :param path_results: Per-path evaluation results to aggregate.
        :type path_results: tuple[PathResult, ...]
        :param total_molecules: Total number of molecules in the input
            dataset, including those without a route. Used as the
            denominator for ``solv_2``.
        :type total_molecules: int

        :return: Aggregated dataset metrics.
        :rtype: DatasetMetrics
        """
        molecules_with_route = len(path_results)
        solved_routes = sum(1 for pr in path_results if pr.is_route_solved)

        all_step_results = [
            sr for pr in path_results for sr in pr.best_variant.step_results
        ]
        total_steps = len(all_step_results)
        passed_steps = sum(1 for sr in all_step_results if sr.score > 0)

        solv_2 = solved_routes / total_molecules if total_molecules > 0 else 0.0

        sum_scores = sum(sr.score for sr in all_step_results)
        mean_chemcensor_score = sum_scores / total_steps if total_steps > 0 else 0.0

        return DatasetMetrics(
            total_molecules=total_molecules,
            molecules_with_route=molecules_with_route,
            solved_routes=solved_routes,
            passed_steps=passed_steps,
            total_steps=total_steps,
            solv_2=solv_2,
            mean_chemcensor_score=mean_chemcensor_score,
        )
