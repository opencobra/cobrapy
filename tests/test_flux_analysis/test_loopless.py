"""Test functionalities of removing loops in model."""

from typing import Callable

import pytest
from optlang.interface import INFEASIBLE, OPTIMAL

from cobra.core import Model
from cobra.flux_analysis.loopless import add_loopless, loopless_solution


@pytest.mark.parametrize("method", ["fastSNP", "potentials", "original"])
def test_loopless_benchmark_before(
    benchmark: Callable, ll_test_model: Model, method: str
) -> None:
    """Benchmark initial condition."""

    def _():
        with ll_test_model:
            add_loopless(ll_test_model, method=method)
            ll_test_model.optimize()

    benchmark(_)


def test_loopless_benchmark_after(benchmark: Callable, ll_test_model: Model) -> None:
    """Benchmark final condition."""
    benchmark(loopless_solution, ll_test_model)


def test_loopless_solution(ll_test_model: Model) -> None:
    """Test loopless_solution()."""
    opt_feasible = ll_test_model.slim_optimize()
    solution_feasible = loopless_solution(ll_test_model)
    ll_test_model.reactions.v3.lower_bound = 1
    opt_infeasible = ll_test_model.slim_optimize()
    solution_infeasible = loopless_solution(ll_test_model)
    assert solution_feasible.fluxes["v3"] == 0.0
    assert solution_feasible.objective_value == pytest.approx(opt_feasible)
    assert solution_infeasible.fluxes["v3"] == 1.0
    assert solution_infeasible.objective_value == pytest.approx(opt_infeasible)


def test_loopless_solution_fluxes(model: Model) -> None:
    """Test fluxes of loopless_solution()."""
    sol = model.optimize()
    ll_solution = loopless_solution(model, fluxes=sol.fluxes)
    assert len(ll_solution.fluxes) == len(model.reactions)
    assert ll_solution.objective_value == pytest.approx(sol.objective_value)


@pytest.mark.parametrize("method", ["fastSNP", "potentials", "original"])
@pytest.mark.parametrize("flux_threshold", [None, 1e-2])
def test_add_loopless(
    ll_test_model: Model, method: str, flux_threshold: float | None
) -> None:
    """Test add_loopless()."""
    add_loopless(ll_test_model, method=method, flux_threshold=flux_threshold)
    feasible_status = ll_test_model.optimize().status
    ll_test_model.reactions.v3.lower_bound = 1
    ll_test_model.slim_optimize()
    infeasible_status = ll_test_model.solver.status
    assert feasible_status == OPTIMAL
    assert infeasible_status == INFEASIBLE
