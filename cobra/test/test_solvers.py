# -*- coding: utf-8 -*-
from __future__ import absolute_import

import pytest

from cobra import solvers
from cobra.core import Metabolite, Model, Reaction

from .conftest import model

try:
    import scipy
except ImportError:
    scipy = None


@pytest.fixture(scope="class", params=list(solvers.solver_dict))
def solver_test(request):
    solver = solvers.solver_dict[request.param]
    old_solution = 0.8739215
    infeasible_model = Model()
    metabolite_1 = Metabolite("met1")
    reaction_1 = Reaction("rxn1")
    reaction_2 = Reaction("rxn2")
    reaction_1.add_metabolites({metabolite_1: 1})
    reaction_2.add_metabolites({metabolite_1: 1})
    reaction_1.lower_bound = 1
    reaction_2.upper_bound = 2
    infeasible_model.add_reactions([reaction_1, reaction_2])
    return solver, old_solution, infeasible_model


class TestCobraSolver:
    def test_attributes(self, solver_test):
        solver, old_solution, infeasible_model = solver_test
        assert hasattr(solver, "create_problem")
        assert hasattr(solver, "solve_problem")
        assert hasattr(solver, "get_status")
        assert hasattr(solver, "get_objective_value")
        assert hasattr(solver, "format_solution")
        assert hasattr(solver, "change_variable_bounds")
        assert hasattr(solver, "change_variable_objective")
        assert hasattr(solver, "solve")
        assert hasattr(solver, "set_parameter")

    def test_creation(self, solver_test, model):
        solver, old_solution, infeasible_model = solver_test
        solver.create_problem(model)

    def test_solve_feasible(self, solver_test, model):
        solver, old_solution, infeasible_model = solver_test
        solution = solver.solve(model)
        assert solution.status == "optimal"
        assert abs(old_solution - solution.f) < 10 ** -4

    def test_solve_minimize(self, solver_test, model):
        solver, old_solution, infeasible_model = solver_test
        solution = solver.solve(model, objective_sense='minimize')
        assert solution.status == "optimal"
        assert abs(solution.f) < 10 ** -4

    def test_low_level_control(self, solver_test):
        solver, old_solution, infeasible_model = solver_test
        lp = solver.create_problem(infeasible_model)
        solver.solve_problem(lp)
        assert solver.get_status(lp) == "infeasible"
        # going to make feasible
        solver.change_variable_bounds(lp, 0, -2., 2.)
        solver.change_variable_bounds(lp, 1, -2., 2.)
        solver.solve_problem(lp)
        # should now be feasible, but obj = 0
        assert solver.get_status(lp) == "optimal"
        assert abs(solver.get_objective_value(lp)) < 10 ** -4
        # should now have obj = 2 (maximize should be the default)
        solver.change_variable_objective(lp, 0, 1.)
        solver.solve_problem(lp)
        assert abs(solver.get_objective_value(lp) - 2) < 10 ** -4
        # should now solve with obj = -2
        solver.solve_problem(lp, objective_sense="minimize")
        assert abs(solver.get_objective_value(lp) + 2) < 10 ** -4
        # should now have obj = 4
        solver.change_variable_objective(lp, 0, 2.)
        solver.solve_problem(lp, objective_sense="maximize")
        assert abs(solver.get_objective_value(lp) - 4) < 10 ** -4
        # make sure the solution looks good still
        solution = solver.format_solution(lp, infeasible_model)
        assert abs(solution.x[0] - 2) < 10 ** -4
        assert abs(solution.x[1] + 2) < 10 ** -4
        assert abs(solution.x_dict["rxn1"] - 2) < 10 ** -4
        assert abs(solution.x_dict["rxn2"] + 2) < 10 ** -4

    def test_set_objective_sense(self, solver_test, model):
        solver, old_solution, infeasible_model = solver_test
        maximize = solver.create_problem(model, objective_sense="maximize")
        minimize = solver.create_problem(model, objective_sense="minimize")
        solver.solve_problem(maximize)
        solver.solve_problem(minimize)
        max_solution = solver.format_solution(maximize, model)
        min_solution = solver.format_solution(minimize, model)
        assert min_solution.status == "optimal"
        assert abs(min_solution.f) < 10 ** -4
        assert abs(old_solution - max_solution.f) < 10 ** -4
        assert max_solution.status == "optimal"
        # if we set minimize at creation, can we override it at solve
        solver.solve_problem(minimize, objective_sense="maximize")
        override_minimize = solver.format_solution(minimize, model)
        assert abs(max_solution.f - override_minimize.f) < 10 ** -4

    def test_solve_mip(self, solver_test):
        solver, old_solution, infeasible_model = solver_test
        if not hasattr(solver, "_SUPPORTS_MILP") or not solver._SUPPORTS_MILP:
            pytest.skip("no milp support")
        cobra_model = Model('MILP_implementation_test')
        constraint = Metabolite("constraint")
        constraint._bound = 2.5
        x = Reaction("x")
        x.lower_bound = 0.
        x.add_metabolites({constraint: 2.5})
        y = Reaction("y")
        y.lower_bound = 0.
        y.add_metabolites({constraint: 1.})
        cobra_model.add_reactions([x, y])
        x.objective_coefficient = 1.
        y.objective_coefficient = 1.
        float_sol = solver.solve(cobra_model)
        # add an integer constraint
        y.variable_kind = "integer"
        int_sol = solver.solve(cobra_model)
        assert abs(float_sol.f - 2.5) < 10 ** -5
        assert abs(float_sol.x_dict["y"] - 2.5) < 10 ** -5
        assert int_sol.status == "optimal"
        assert abs(int_sol.f - 2.2) < 10 ** -3
        assert abs(int_sol.x_dict["y"] - 2.0) < 10 ** -3

    def test_solve_infeasible(self, solver_test):
        solver, old_solution, infeasible_model = solver_test
        solution = solver.solve(infeasible_model)
        assert solution.status == "infeasible"

    def test_independent_creation(self, solver_test, model):
        solver, old_solution, infeasible_model = solver_test
        feasible_lp = solver.create_problem(model)
        infeasible_lp = solver.create_problem(infeasible_model)
        solver.solve_problem(feasible_lp)
        solver.solve_problem(infeasible_lp)
        feasible_solution = solver.format_solution(feasible_lp, model)
        infeasible_solution = solver.format_solution(infeasible_lp,
                                                     infeasible_model)
        assert feasible_solution.status == "optimal"
        assert abs(old_solution - feasible_solution.f) < 10 ** -4
        assert infeasible_solution.status == "infeasible"

    def test_change_coefficient(self, solver_test):
        solver, old_solution, infeasible_model = solver_test
        c = Metabolite("c")
        c._bound = 6
        x = Reaction("x")
        x.lower_bound = 1.
        y = Reaction("y")
        y.lower_bound = 0.
        x.add_metabolites({c: 1})
        z = Reaction("z")
        z.add_metabolites({c: 1})
        m = Model("test_model")
        m.add_reactions([x, y, z])
        z.objective_coefficient = 1
        # change an existing coefficient
        lp = solver.create_problem(m)
        solver.solve_problem(lp)
        sol1 = solver.format_solution(lp, m)
        assert sol1.status == "optimal"
        solver.change_coefficient(lp, 0, 0, 2)
        solver.solve_problem(lp)
        sol2 = solver.format_solution(lp, m)
        assert sol2.status == "optimal"
        assert abs(sol1.f - 5.0) < 10 ** -3
        assert abs(sol2.f - 4.0) < 10 ** -3
        # change a new coefficient
        z.objective_coefficient = 0.
        y.objective_coefficient = 1.
        lp = solver.create_problem(m)
        solver.change_coefficient(lp, 0, 1, 2)
        solver.solve_problem(lp)
        solution = solver.format_solution(lp, m)
        assert solution.status == "optimal"
        assert abs(solution.x_dict["y"] - 2.5) < 10 ** -3

    def test_inequality(self, solver_test):
        solver, old_solution, infeasible_model = solver_test
        # The space enclosed by the constraints is a 2D triangle with
        # vertexes as (3, 0), (1, 2), and (0, 1)
        # c1 encodes y - x > 1 ==> y > x - 1
        # c2 encodes y + x < 3 ==> y < 3 - x
        c1 = Metabolite("c1")
        c2 = Metabolite("c2")
        x = Reaction("x")
        x.lower_bound = 0
        y = Reaction("y")
        y.lower_bound = 0
        x.add_metabolites({c1: -1, c2: 1})
        y.add_metabolites({c1: 1, c2: 1})
        c1._bound = 1
        c1._constraint_sense = "G"
        c2._bound = 3
        c2._constraint_sense = "L"
        m = Model()
        m.add_reactions([x, y])
        # test that optimal values are at the vertices
        m.objective = "x"
        assert abs(solver.solve(m).f - 1.0) < 10 ** -3
        assert abs(solver.solve(m).x_dict["y"] - 2.0) < 10 ** -3
        m.objective = "y"
        assert abs(solver.solve(m).f - 3.0) < 10 ** -3
        assert abs(
            solver.solve(m, objective_sense="minimize").f - 1.0) < 10 ** -3

    @pytest.mark.skipif(scipy is None,
                        reason="scipy required for quadratic objectives")
    def test_quadratic(self, solver_test):
        solver, old_solution, infeasible_model = solver_test
        if not hasattr(solver, "set_quadratic_objective"):
            pytest.skip("no qp support")
        c = Metabolite("c")
        c._bound = 2
        x = Reaction("x")
        x.lower_bound = 0.
        y = Reaction("y")
        y.lower_bound = 0.
        x.add_metabolites({c: 1})
        y.add_metabolites({c: 1})
        m = Model()
        m.add_reactions([x, y])
        x.objective_coefficient = -0.5
        y.objective_coefficient = -0.5
        lp = solver.create_problem(m)
        quadratic_obj = scipy.sparse.eye(2) * 2
        solver.set_quadratic_objective(lp, quadratic_obj)
        solver.solve_problem(lp, objective_sense="minimize")
        solution = solver.format_solution(lp, m)
        assert solution.status == "optimal"
        # Respecting linear objectives also makes the objective value 1.
        assert abs(solution.f - 1.) < 10 ** -3
        assert abs(solution.x_dict["y"] - 1.) < 10 ** -3
        assert abs(solution.x_dict["y"] - 1.) < 10 ** -3
        # When the linear objectives are removed the objective value is 2.
        solver.change_variable_objective(lp, 0, 0.)
        solver.change_variable_objective(lp, 1, 0.)
        solver.solve_problem(lp, objective_sense="minimize")
        solution = solver.format_solution(lp, m)
        assert solution.status == "optimal"
        assert abs(solution.f - 2.) < 10 ** -3
        # test quadratic from solve function
        solution = solver.solve(m, quadratic_component=quadratic_obj,
                                objective_sense="minimize")
        assert solution.status == "optimal"
        assert abs(solution.f - 1.) < 10 ** -3
        c._bound = 6
        z = Reaction("z")
        x.objective_coefficient = 0.
        y.objective_coefficient = 0.
        z.lower_bound = 0.
        z.add_metabolites({c: 1})
        m.add_reaction(z)
        solution = solver.solve(m, quadratic_component=scipy.sparse.eye(3),
                                objective_sense="minimize")
        # should be 12 not 24 because 1/2 (V^T Q V)
        assert solution.status == "optimal"
        assert abs(solution.f - 6) < 10 ** -3
        assert abs(solution.x_dict["x"] - 2) < 10 ** -6
        assert abs(solution.x_dict["y"] - 2) < 10 ** -6
        assert abs(solution.x_dict["z"] - 2) < 10 ** -6
