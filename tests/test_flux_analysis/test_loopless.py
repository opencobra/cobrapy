"""Test functionalities of removing loops in model."""

from typing import Callable

import pytest
from optlang.interface import INFEASIBLE, OPTIMAL

from cobra.core import Metabolite, Model, Reaction
from cobra.flux_analysis.loopless import add_loopless, loopless_solution
from cobra.util import solver as sutil


def construct_ll_test_model() -> Model:
    """Construct test model."""
    test_model = Model()
    test_model.add_metabolites(Metabolite("A"))
    test_model.add_metabolites(Metabolite("B"))
    test_model.add_metabolites(Metabolite("C"))
    EX_A = Reaction("EX_A")
    EX_A.add_metabolites({test_model.metabolites.A: 1})
    DM_C = Reaction("DM_C")
    DM_C.add_metabolites({test_model.metabolites.C: -1})
    v1 = Reaction("v1")
    v1.add_metabolites({test_model.metabolites.A: -1, test_model.metabolites.B: 1})
    v2 = Reaction("v2")
    v2.add_metabolites({test_model.metabolites.B: -1, test_model.metabolites.C: 1})
    v3 = Reaction("v3")
    v3.add_metabolites({test_model.metabolites.C: -1, test_model.metabolites.A: 1})
    test_model.add_reactions([EX_A, DM_C, v1, v2, v3])
    DM_C.objective_coefficient = 1
    return test_model


@pytest.fixture(
    scope="function",
    params=[s for s in ["glpk", "cplex", "gurobi"] if s in sutil.solvers],
)
def ll_test_model(request: pytest.FixtureRequest) -> Model:
    """Return test model set with different solvers."""
    test_model = construct_ll_test_model()
    test_model.solver = request.param
    return test_model


def test_loopless_benchmark_before(benchmark: Callable) -> None:
    """Benchmark initial condition."""
    test_model = construct_ll_test_model()

    def _():
        with test_model:
            add_loopless(test_model)
            test_model.optimize()

    benchmark(_)


def test_loopless_benchmark_after(benchmark: Callable) -> None:
    """Benchmark final condition."""
    test_model = construct_ll_test_model()
    benchmark(loopless_solution, test_model)


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


def test_add_loopless(ll_test_model: Model) -> None:
    """Test add_loopless()."""
    add_loopless(ll_test_model)
    feasible_status = ll_test_model.optimize().status
    ll_test_model.reactions.v3.lower_bound = 1
    ll_test_model.slim_optimize()
    infeasible_status = ll_test_model.solver.status
    assert feasible_status == OPTIMAL
    assert infeasible_status == INFEASIBLE
