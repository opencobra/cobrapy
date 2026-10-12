"""Test functionalities of pFBA."""

import warnings
from typing import Callable, List

import pytest

from cobra.core import Model
from cobra.exceptions import Infeasible, OptimizationError
from cobra.flux_analysis import parsimonious
from cobra.flux_analysis.parsimonious import add_pfba, pfba


def test_pfba_benchmark(
    large_model: Model, benchmark: Callable, all_solvers: List[str]
) -> None:
    """Benchmark pFBA functionality."""
    large_model.solver = all_solvers
    benchmark(pfba, large_model)


def test_pfba(model: Model, all_solvers: List[str]) -> None:
    """Test pFBA functionality."""
    model.solver = all_solvers
    with model:
        add_pfba(model)
        with pytest.raises(ValueError):
            add_pfba(model)

    expression = model.objective.expression
    n_constraints = len(model.constraints)
    solution = pfba(model)
    assert solution.status == "optimal"
    assert solution.fluxes["Biomass_Ecoli_core"] == pytest.approx(
        0.8739, abs=1e-4, rel=0.0
    )
    assert solution.fluxes.abs().sum() == pytest.approx(518.4221, abs=1e-4, rel=0.0)

    # test changes to model reverted
    assert expression == model.objective.expression
    assert len(model.constraints) == n_constraints

    # needed?
    # Test desired_objective_value
    # desired_objective = 0.8
    # pfba(model, solver=solver,
    #                       desired_objective_value=desired_objective)
    # abs_x = [abs(i) for i in model.solution.x]
    # assert model.solution.status == "optimal"
    # assert abs(model.solution.f - desired_objective) < 0.001
    # assert abs(sum(abs_x) - 476.1594) < 0.001

    # TODO: parametrize fraction (DRY it up)
    # Test fraction_of_optimum
    solution = pfba(model, fraction_of_optimum=0.95)
    assert solution.status == "optimal"
    assert solution.fluxes["Biomass_Ecoli_core"] == pytest.approx(
        0.95 * 0.8739, abs=1e-4, rel=0.0
    )
    abs_x = [abs(i) for i in solution.fluxes.values]
    assert sum(abs_x) == pytest.approx(493.4400, abs=1e-4, rel=0.0)

    # Infeasible solution
    model.reactions.ATPM.lower_bound = 500
    with warnings.catch_warnings():
        warnings.simplefilter("error", UserWarning)
        with pytest.raises((UserWarning, Infeasible, ValueError)):
            pfba(model)


@pytest.mark.parametrize("raise_error", [True, False])
def test_pfba_forwards_raise_error(model: Model, mocker, raise_error: bool) -> None:
    """Test that `raise_error` is forwarded to `get_solution` like in `optimize`."""
    get_solution = mocker.spy(parsimonious, "get_solution")
    solution = pfba(model, raise_error=raise_error)
    assert solution.status == "optimal"
    assert get_solution.call_args.kwargs["raise_error"] is raise_error


def test_pfba_raise_error_infeasible(model: Model) -> None:
    """Test that an infeasible problem raises when `raise_error` is True."""
    with model:
        model.reactions.ATPM.lower_bound = 500
        with pytest.raises(OptimizationError):
            pfba(model, raise_error=True)
