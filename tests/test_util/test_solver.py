"""Test functions of solver.py."""

import logging
from typing import TYPE_CHECKING, Optional

import numpy as np
import pytest

from cobra.exceptions import OptimizationError
from cobra.util import solver as su


if TYPE_CHECKING:
    from cobra import Model

stable_optlang = ["glpk", "cplex", "gurobi"]
optlang_solvers = [f"optlang-{s}" for s in stable_optlang if s in su.solvers]


def test_solver_list() -> None:
    """Expect that at least the GLPK solver is found."""
    assert len(su.solvers) >= 1
    assert "glpk" in su.solvers


def test_interface_str() -> None:
    """Test the string representation of solver interfaces."""
    assert su.interface_to_str("nonsense") == "nonsense"
    assert su.interface_to_str("optlang.glpk_interface") == "glpk"
    assert su.interface_to_str("optlang-cplex") == "cplex"


def test_solver_name() -> None:
    """Test that the default LP solver name is GLPK."""
    assert su.get_solver_name() == "glpk"


def test_choose_solver(model: "Model") -> Optional[su.SolverNotFound]:
    """Test that solver switching is working."""
    so = su.choose_solver(model, "glpk")
    assert su.interface_to_str(so) == "glpk"

    if any(s in su.solvers for s in su.qp_solvers):
        qp_choice = su.choose_solver(model, qp=True)
        assert su.interface_to_str(qp_choice) in su.qp_solvers
    else:
        with pytest.raises(su.SolverNotFound):
            su.choose_solver(model, qp=True)


def test_linear_reaction_coefficients(model: "Model") -> None:
    """Test that linear coefficients are identifiable in objective."""
    coefficients = su.linear_reaction_coefficients(model)
    assert coefficients == {model.reactions.Biomass_Ecoli_core: 1}


def test_fail_non_linear_reaction_coefficients(model: "Model") -> None:
    """Test failure of non-linear coefficient identification in reaction."""
    model.solver = "optlang-glpk"

    with pytest.raises(ValueError) as error:
        model.objective = model.problem.Objective(
            model.reactions.ATPM.flux_expression**2
        )
        coefficients = su.linear_reaction_coefficients(model)
        assert coefficients == {}

    assert "GLPK only supports linear objectives." in str(error.value)


def test_add_remove(model: "Model") -> None:
    """Test addition and removal of variables and constraints."""
    v = model.variables
    new_var = model.problem.Variable("test_var", lb=-10, ub=-10)
    new_constraint = model.problem.Constraint(
        v.PGK - new_var, name="test_constraint", lb=0
    )

    su.add_cons_vars_to_problem(model, [new_var, new_constraint])
    assert "test_var" in model.variables.keys()
    assert "test_constraint" in model.constraints.keys()

    su.remove_cons_vars_from_problem(model, [new_var, new_constraint])
    assert "test_var" not in model.variables.keys()
    assert "test_constraint" not in model.constraints.keys()


def test_add_remove_in_context(model: "Model") -> None:
    """Test addition and removal of variables and constraints within context."""
    v = model.variables
    new_var = model.problem.Variable("test_var", lb=-10, ub=-10)

    with model:
        su.add_cons_vars_to_problem(model, [new_var])
        su.remove_cons_vars_from_problem(model, [v.PGM])
        assert "test_var" in model.variables.keys()
        assert "PGM" not in model.variables.keys()

    assert "test_var" not in model.variables.keys()
    assert "PGM" in model.variables.keys()


def test_absolute_expression(model: "Model") -> None:
    """Test addition of an absolute expression."""
    v = model.variables
    with model:
        parts = su.add_absolute_expression(model, 2 * v.PGM, name="test", ub=100)
        assert len(parts) == 3
        assert "test" in model.variables.keys()
        assert "abs_pos_test" in model.constraints.keys()
        assert "abs_neg_test" in model.constraints.keys()
    assert "test" not in model.variables.keys()
    assert "abs_pos_test" not in model.constraints.keys()
    assert "abs_neg_test" not in model.constraints.keys()


@pytest.mark.parametrize("solver", optlang_solvers)
def test_fix_objective_as_constraint(solver: str, model: "Model") -> None:
    """Test fixing present objective as a constraint."""
    model.solver = solver
    opt = model.slim_optimize()
    with model as m:
        su.fix_objective_as_constraint(model, 1.0, name="fixed")
        assert (m.constraints.fixed.expression - m.objective.expression).simplify() == 0
        assert m.constraints.fixed.lb == pytest.approx(opt)
    assert "fixed" not in m.constraints
    su.fix_objective_as_constraint(model, name="fixed")
    assert (
        model.constraints.fixed.expression - model.objective.expression
    ).simplify() == 0
    assert m.constraints.fixed.lb == pytest.approx(opt)
    assert "fixed" in model.constraints


@pytest.mark.parametrize("solver", optlang_solvers)
def test_fix_objective_as_constraint_minimize(model: "Model", solver: str) -> None:
    """Test fixing present objective as a constraint but as a minimization."""
    model.solver = solver
    model.reactions.Biomass_Ecoli_core.bounds = (0.1, 0.1)
    minimize_glucose = model.problem.Objective(
        model.reactions.EX_glc__D_e.flux_expression, direction="min"
    )
    su.set_objective(model, minimize_glucose)
    su.fix_objective_as_constraint(model)
    fx_name = f"fixed_objective_{model.objective.name}"
    constr = model.constraints
    # Ensure that a solution exists on non-GLPK solvers.
    model.slim_optimize()
    assert (constr[fx_name].lb, constr[fx_name].ub) == (
        None,
        model.solver.objective.value,
    )


@pytest.mark.parametrize("solver", optlang_solvers)
def test_add_lp_feasibility(model: "Model", solver: str) -> None:
    """Test functionality to ensure LP feasibility."""
    model.solver = solver

    with model:
        with model:
            su.add_lp_feasibility(model)
            assert "s_plus_succoa_c" in model.variables
            assert np.isclose(model.slim_optimize(), 0.0)

        model.reactions.EX_glc__D_e.lower_bound = 1
        assert np.isnan(model.slim_optimize())
        assert "s_plus_succoa_c" not in model.variables

        su.add_lp_feasibility(model)
        assert np.isclose(model.slim_optimize(), 1.3872307692307695)


@pytest.mark.parametrize("solver", optlang_solvers)
def test_add_lexicographic_constraints(model: "Model", solver: str) -> None:
    """Test addition of lexicographic constraints."""
    model.solver = solver

    rxns = ["Biomass_Ecoli_core", "EX_glc__D_e", "EX_o2_e"]

    with model:
        out = su.add_lexicographic_constraints(model, rxns, ["max", "min", "max"])
        print(model.reactions.Biomass_Ecoli_core.bounds)
        assert np.isclose(model.constraints[-3].lb, out[rxns[0]])
        assert np.isclose(model.constraints[-2].ub, out[rxns[1]])
        assert np.isclose(model.constraints[-1].lb, out[rxns[2]])

    with model:
        su.add_lexicographic_constraints(model, rxns, "max")

    with model:
        su.add_lexicographic_constraints(model, rxns)


def test_time_limit(large_model: "Model") -> None:
    """Test time limit while optimizing a model."""
    if su.interface_to_str(large_model.problem) != "glpk":
        pytest.skip("requires GLPK")

    # It is done like this since optlang accepts inputs in seconds
    # whereas GLPK accepts milliseconds
    large_model.solver.configuration._smcp.tm_lim = 1
    with pytest.warns(UserWarning):
        sol = large_model.optimize()
    assert sol.fluxes is not None

    with pytest.raises(OptimizationError):
        sol = large_model.optimize(raise_error=True)


@pytest.mark.parametrize(
    "solver", [s for s in su.solvers if s in ["osqp", "coinor_cbc"]]
)
def test_specialized_solver_warning(solver, caplog):
    """Test the warning for specialized solvers."""
    with caplog.at_level(logging.WARNING):
        su.check_solver(solver)
    assert "are specialized solvers for" in caplog.text
