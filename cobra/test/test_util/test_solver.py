# -*- coding: utf-8 -*-

"""Test functions of solver.py"""

from __future__ import absolute_import

import pytest

import cobra.util.solver as su
from cobra.exceptions import OptimizationError

stable_optlang = ["glpk", "cplex", "gurobi"]
optlang_solvers = ["optlang-" + s for s in stable_optlang if s in su.solvers]


def test_solver_list():
    assert len(su.solvers) >= 1
    assert "glpk" in su.solvers


def test_interface_str():
    assert su.interface_to_str("nonsense") == "nonsense"
    assert su.interface_to_str("optlang.glpk_interface") == "glpk"
    assert su.interface_to_str("optlang-cplex") == "cplex"


def test_solver_name():
    assert su.get_solver_name() == "glpk"


def test_choose_solver(model):
    so = su.choose_solver(model)
    assert su.interface_to_str(so) == "glpk"
    so = su.choose_solver(model, "glpk")
    assert su.interface_to_str(so) == "glpk"

    if any(s in su.solvers for s in su.qp_solvers):
        qp_choice = su.choose_solver(model, qp=True)
        assert su.interface_to_str(qp_choice) in su.qp_solvers
    else:
        with pytest.raises(su.SolverNotFound):
            su.choose_solver(model, qp=True)


def test_linear_reaction_coefficients(model):
    coefficients = su.linear_reaction_coefficients(model)
    assert coefficients == {model.reactions.Biomass_Ecoli_core: 1}


@pytest.mark.parametrize("solver", optlang_solvers)
def test_fail_non_linear_reaction_coefficients(model, solver):
    model.solver = solver
    try:
        model.objective = model.problem.Objective(
            model.reactions.ATPM.flux_expression ** 2
        )
    except ValueError:
        pass
    else:
        coefficients = su.linear_reaction_coefficients(model)
        assert coefficients == {}
        with pytest.raises(ValueError):
            model.reactions.ACALD.objective_coefficient = 1


def test_add_remove(model):
    v = model.variables
    new_var = model.problem.Variable("test_var", lb=-10, ub=-10)
    new_constraint = model.problem.Constraint(
        v.PGK - new_var, name="test_constraint", lb=0)

    su.add_cons_vars_to_problem(model, [new_var, new_constraint])
    assert "test_var" in model.variables.keys()
    assert "test_constraint" in model.constraints.keys()

    su.remove_cons_vars_from_problem(model, [new_var, new_constraint])
    assert "test_var" not in model.variables.keys()
    assert "test_constraint" not in model.constraints.keys()


def test_add_remove_in_context(model):
    v = model.variables
    new_var = model.problem.Variable("test_var", lb=-10, ub=-10)

    with model:
        su.add_cons_vars_to_problem(model, [new_var])
        su.remove_cons_vars_from_problem(model, [v.PGM])
        assert "test_var" in model.variables.keys()
        assert "PGM" not in model.variables.keys()

    assert "test_var" not in model.variables.keys()
    assert "PGM" in model.variables.keys()


def test_absolute_expression(model):
    v = model.variables
    with model:
        parts = su.add_absolute_expression(
            model, 2 * v.PGM, name="test", ub=100)
        assert len(parts) == 3
        assert "test" in model.variables.keys()
        assert "abs_pos_test" in model.constraints.keys()
        assert "abs_neg_test" in model.constraints.keys()
    assert "test" not in model.variables.keys()
    assert "abs_pos_test" not in model.constraints.keys()
    assert "abs_neg_test" not in model.constraints.keys()


@pytest.mark.parametrize("solver", optlang_solvers)
def test_fix_objective_as_constraint(solver, model):
    model.solver = solver
    with model as m:
        su.fix_objective_as_constraint(model, 1.0)
        constraint_name = m.constraints[-1]
        assert abs(m.constraints[-1].expression -
                   m.objective.expression) < 1e-6
    assert constraint_name not in m.constraints
    su.fix_objective_as_constraint(model)
    constraint_name = model.constraints[-1]
    assert abs(model.constraints[-1].expression -
               model.objective.expression) < 1e-6
    assert constraint_name in model.constraints


@pytest.mark.parametrize("solver", optlang_solvers)
def test_fix_objective_as_constraint_minimize(model, solver):
    model.solver = solver
    model.reactions.Biomass_Ecoli_core.bounds = (0.1, 0.1)
    minimize_glucose = model.problem.Objective(
        model.reactions.EX_glc__D_e.flux_expression,
        direction='min')
    su.set_objective(model, minimize_glucose)
    su.fix_objective_as_constraint(model)
    fx_name = 'fixed_objective_{}'.format(model.objective.name)
    constr = model.constraints
    # Ensure that a solution exists on non-GLPK solvers.
    model.slim_optimize()
    assert (constr[fx_name].lb, constr[fx_name].ub) == (
        None, model.solver.objective.value)


def test_time_limit(large_model):
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
