import pytest

import cobra.util.solver as su
from cobra.test.conftest import model
import optlang

stable_optlang = ["glpk", "cplex", "gurobi"]
optlang_solvers = ["optlang-" + s for s in stable_optlang if s in su.solvers]


class TestHelpers:
    def test_solver_list(self):
        assert len(su.solvers) >= 1
        assert "glpk" in su.solvers

    def test_interface_str(self):
        assert su.interface_to_str("nonsense") == "nonsense"
        assert su.interface_to_str("optlang.glpk_interface") == "glpk"
        assert su.interface_to_str("optlang-cplex") == "cplex"

    def test_solver_name(self):
        assert su.get_solver_name() == "glpk"

    def test_choose_solver(self, model):
        legacy, so = su.choose_solver(model)
        assert not legacy
        assert su.interface_to_str(so) == "glpk"
        legacy, so = su.choose_solver(model, "optlang-glpk")
        assert not legacy
        assert su.interface_to_str(so) == "glpk"
        assert su.choose_solver(model, "cglpk")[0]

        if any(s in su.solvers for s in su.qp_solvers):
            qp_choice = su.interface_to_str(su.choose_solver(model, qp=True))
            assert qp_choice in su.qp_solvers
        else:
            with pytest.raises(su.SolverNotFound):
                su.choose_solver(model, qp=True)


class TestObjectiveHelpers:
    def test_linear_reaction_coefficients(self, model):
        coefficients = su.linear_reaction_coefficients(model)
        assert coefficients == {model.reactions.Biomass_Ecoli_core: 1}

    @pytest.mark.parametrize("solver", optlang_solvers)
    def test_fail_non_linear_reaction_coefficients(self, model, solver):
        model.solver = solver
        try:
            model.objective = model.solver.interface.Objective(
                model.reactions.ATPM.flux_expression ** 2
            )
        except ValueError:
            pass
        else:
            coefficients = su.linear_reaction_coefficients(model)
            assert coefficients == {}
            with pytest.raises(ValueError):
                model.reactions.ACALD.objective_coefficient = 1


class TestSolverMods:
    def test_add_remove(self, model):
        v = model.solver.variables
        new_var = model.solver.interface.Variable("test_var", lb=-10, ub=-10)
        new_constraint = model.solver.interface.Constraint(
            v.PGK - new_var, name="test_constraint", lb=0)

        su.add_to_solver(model, [new_var, new_constraint])
        assert "test_var" in model.solver.variables.keys()
        assert "test_constraint" in model.solver.constraints.keys()

        su.remove_from_solver(model, [new_var, new_constraint])
        assert "test_var" not in model.solver.variables.keys()
        assert "test_constraint" not in model.solver.constraints.keys()

    def test_add_remove_in_context(self, model):
        v = model.solver.variables
        new_var = model.solver.interface.Variable("test_var", lb=-10, ub=-10)

        with model:
            su.add_to_solver(model, [new_var])
            su.remove_from_solver(model, [v.PGM])
            assert "test_var" in model.solver.variables.keys()
            assert "PGM" not in model.solver.variables.keys()

        assert "test_var" not in model.solver.variables.keys()
        assert "PGM" in model.solver.variables.keys()

    def test_absolute_expression(self, model):
        v = model.solver.variables
        with model:
            su.add_absolute_expression(model, 2 * v.PGM, name="test", ub=100)
            assert "test" in model.solver.variables.keys()
            assert "abs_pos_test" in model.solver.constraints.keys()
            assert "abs_neg_test" in model.solver.constraints.keys()
        assert "test" not in model.solver.variables.keys()
        assert "abs_pos_test" not in model.solver.constraints.keys()
        assert "abs_neg_test" not in model.solver.constraints.keys()

    @pytest.mark.parametrize("solver", optlang_solvers)
    def test_fix_objective_as_constraint(self, solver, model):
        model.solver = solver
        with model as m:
            su.fix_objective_as_constraint(model, 1.0)
            constraint_name = m.solver.constraints[-1]
            assert abs(m.solver.constraints[-1].expression -
                       m.objective.expression) < 1e-6
        assert constraint_name not in m.solver.constraints
        su.fix_objective_as_constraint(model)
        constraint_name = model.solver.constraints[-1]
        assert abs(model.solver.constraints[-1].expression -
                   model.objective.expression) < 1e-6
        assert constraint_name in model.solver.constraints

    @pytest.mark.parametrize("solver", optlang_solvers)
    def test_fix_objective_as_constraint_minimize(self, solver, model):
        model.reactions.Biomass_Ecoli_core.bounds = (0.1, 0.1)
        minimize_glucose = model.solver.interface.Objective(
            model.reactions.EX_glc__D_e.flux_expression,
            direction='min')
        su.set_objective(model, minimize_glucose)
        su.fix_objective_as_constraint(model)
        fx_name = 'Fixed_objective_{}'.format(model.objective.name)
        constr = model.solver.constraints
        assert (constr[fx_name].lb, constr[fx_name].ub) == (
            None, model.solution.objective_value)
