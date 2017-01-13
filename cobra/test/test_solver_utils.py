import cobra.util.solver as su
from cobra.test.conftest import model


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
