from unittest import TestCase, TestLoader, TextTestRunner, skipIf
import sys
# deal with absolute imports by adding the appropriate directory to the path
if __name__ == "__main__":
    sys.path.insert(0, "../..")
    from cobra.test import create_test_model
    from cobra import Model, Reaction, Metabolite
    from cobra import solvers
    sys.path.pop(0)  # remove the added directory to the path
else:
    from . import create_test_model
    from .. import Model, Reaction, Metabolite
    from .. import solvers

try:
    import scipy
except:
    scipy = None


class TestCobraSolver(object):
    def setUp(self):
        self.solver = solvers.solver_dict[self.solver_name]
        self.model = create_test_model("textbook")
        self.old_solution = 0.8739215
        self.infeasible_model = Model()
        metabolite_1 = Metabolite("met1")
        reaction_1 = Reaction("rxn1")
        reaction_2 = Reaction("rxn2")
        reaction_1.add_metabolites({metabolite_1: 1})
        reaction_2.add_metabolites({metabolite_1: 1})
        reaction_1.lower_bound = 1
        reaction_2.upper_bound = 2
        self.infeasible_model.add_reactions([reaction_1, reaction_2])

    def test_attributes(self):
        solver = self.solver
        self.assertTrue(hasattr(solver, "create_problem"))
        self.assertTrue(hasattr(solver, "solve_problem"))
        self.assertTrue(hasattr(solver, "get_status"))
        self.assertTrue(hasattr(solver, "get_objective_value"))
        self.assertTrue(hasattr(solver, "format_solution"))
        self.assertTrue(hasattr(solver, "change_variable_bounds"))
        self.assertTrue(hasattr(solver, "change_variable_objective"))
        self.assertTrue(hasattr(solver, "solve"))
        self.assertTrue(hasattr(solver, "set_parameter"))
        # self.assertTrue(hasattr(solver, "update_problem"))

    def test_creation(self):
        solver = self.solver
        solver.create_problem(self.model)

    def test_solve_feasible(self):
        solver = self.solver
        solution = solver.solve(self.model)
        self.assertEqual(solution.status, "optimal")
        self.assertAlmostEqual(self.old_solution,
                               solution.f, places=4)

    def test_solve_minimize(self):
        solver = self.solver
        solution = solver.solve(self.model, objective_sense='minimize')
        self.assertEqual(solution.status, "optimal")
        self.assertAlmostEqual(0, solution.f, places=4)

    def test_low_level_control(self):
        solver = self.solver
        lp = solver.create_problem(self.infeasible_model)
        solver.solve_problem(lp)
        self.assertEqual(solver.get_status(lp), "infeasible")
        # going to make feasible
        solver.change_variable_bounds(lp, 0, -2., 2.)
        solver.change_variable_bounds(lp, 1, -2., 2.)
        solver.solve_problem(lp)
        # should now be feasible, but obj = 0
        self.assertEqual(solver.get_status(lp), "optimal")
        self.assertAlmostEqual(solver.get_objective_value(lp), 0, places=4)
        # should now have obj = 2 (maximize should be the default)
        solver.change_variable_objective(lp, 0, 1.)
        solver.solve_problem(lp)
        self.assertAlmostEqual(solver.get_objective_value(lp), 2, places=4)
        # should now solve with obj = -2
        solver.solve_problem(lp, objective_sense="minimize")
        self.assertAlmostEqual(solver.get_objective_value(lp), -2, places=4)
        # should now have obj = 4
        solver.change_variable_objective(lp, 0, 2.)
        solver.solve_problem(lp, objective_sense="maximize")
        self.assertAlmostEqual(solver.get_objective_value(lp), 4, places=4)
        # make sure the solution looks good still
        solution = solver.format_solution(lp, self.infeasible_model)
        self.assertAlmostEqual(solution.x[0], 2, places=4)
        self.assertAlmostEqual(solution.x[1], -2, places=4)
        self.assertAlmostEqual(solution.x_dict["rxn1"], 2, places=4)
        self.assertAlmostEqual(solution.x_dict["rxn2"], -2, places=4)

    def test_set_objective_sense(self):
        solver = self.solver
        maximize = solver.create_problem(self.model,
                                         objective_sense="maximize")
        minimize = solver.create_problem(self.model,
                                         objective_sense="minimize")
        solver.solve_problem(maximize)
        solver.solve_problem(minimize)
        max_solution = solver.format_solution(maximize, self.model)
        min_solution = solver.format_solution(minimize, self.model)
        self.assertEqual(min_solution.status, "optimal")
        self.assertAlmostEqual(0, min_solution.f, places=4)
        self.assertAlmostEqual(self.old_solution, max_solution.f, places=4)
        self.assertEqual(max_solution.status, "optimal")
        # if we set minimize at creation, can we override it at solve
        solver.solve_problem(minimize, objective_sense="maximize")
        override_minimize = solver.format_solution(minimize, self.model)
        self.assertAlmostEqual(max_solution.f, override_minimize.f, places=4)

    def test_solve_mip(self):
        solver = self.solver
        if not hasattr(solver, "_SUPPORTS_MILP") or not solver._SUPPORTS_MILP:
            self.skipTest("no milp support")
        cobra_model = Model('MILP_implementation_test')
        constraint = Metabolite("constraint")
        constraint._bound = 2.5
        x = Reaction("x")
        x.lower_bound = 0.
        x.objective_coefficient = 1.
        x.add_metabolites({constraint: 2.5})
        y = Reaction("y")
        y.lower_bound = 0.
        y.objective_coefficient = 1.
        y.add_metabolites({constraint: 1.})
        cobra_model.add_reactions([x, y])
        float_sol = solver.solve(cobra_model)
        # add an integer constraint
        y.variable_kind = "integer"
        int_sol = solver.solve(cobra_model)
        self.assertAlmostEqual(float_sol.f, 2.5)
        self.assertAlmostEqual(float_sol.x_dict["y"], 2.5)
        self.assertEqual(int_sol.status, "optimal")
        self.assertAlmostEqual(int_sol.f, 2.2)
        self.assertAlmostEqual(int_sol.x_dict["y"], 2.0)

    def test_solve_infeasible(self):
        solver = self.solver
        solution = solver.solve(self.infeasible_model)
        self.assertEqual(solution.status, "infeasible")

    def test_independent_creation(self):
        solver = self.solver
        feasible_lp = solver.create_problem(self.model)
        infeasible_lp = solver.create_problem(self.infeasible_model)
        solver.solve_problem(feasible_lp)
        solver.solve_problem(infeasible_lp)
        feasible_solution = solver.format_solution(feasible_lp, self.model)
        infeasible_solution = solver.format_solution(infeasible_lp,
                                                     self.infeasible_model)
        self.assertEqual(feasible_solution.status, "optimal")
        self.assertAlmostEqual(self.old_solution,
                               feasible_solution.f, places=4)
        self.assertEqual(infeasible_solution.status, "infeasible")

    def test_change_coefficient(self):
        solver = self.solver
        c = Metabolite("c")
        c._bound = 6
        x = Reaction("x")
        x.lower_bound = 1.
        y = Reaction("y")
        y.lower_bound = 0.
        x.add_metabolites({c: 1})
        z = Reaction("z")
        z.add_metabolites({c: 1})
        z.objective_coefficient = 1
        m = Model("test_model")
        m.add_reactions([x, y, z])
        # change an existing coefficient
        lp = solver.create_problem(m)
        solver.solve_problem(lp)
        sol1 = solver.format_solution(lp, m)
        self.assertEqual(sol1.status, "optimal")
        solver.change_coefficient(lp, 0, 0, 2)
        solver.solve_problem(lp)
        sol2 = solver.format_solution(lp, m)
        self.assertEqual(sol2.status, "optimal")
        self.assertAlmostEqual(sol1.f, 5.0)
        self.assertAlmostEqual(sol2.f, 4.0)
        # change a new coefficient
        z.objective_coefficient = 0.
        y.objective_coefficient = 1.
        lp = solver.create_problem(m)
        solver.change_coefficient(lp, 0, 1, 2)
        solver.solve_problem(lp)
        solution = solver.format_solution(lp, m)
        self.assertEqual(solution.status, "optimal")
        self.assertAlmostEqual(solution.x_dict["y"], 2.5)

    def test_inequality(self):
        # The space enclosed by the constraints is a 2D triangle with
        # vertexes as (3, 0), (1, 2), and (0, 1)
        solver = self.solver
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
        self.assertAlmostEqual(solver.solve(m).f, 1.0)
        self.assertAlmostEqual(solver.solve(m).x_dict["y"], 2.0)
        m.objective = "y"
        self.assertAlmostEqual(solver.solve(m).f, 3.0)
        self.assertAlmostEqual(solver.solve(m, objective_sense="minimize").f,
                               1.0)

    @skipIf(scipy is None, "scipy required for quadratic objectives")
    def test_quadratic(self):
        solver = self.solver
        if not hasattr(solver, "set_quadratic_objective"):
            self.skipTest("no qp support")
        c = Metabolite("c")
        c._bound = 2
        x = Reaction("x")
        x.objective_coefficient = -0.5
        x.lower_bound = 0.
        y = Reaction("y")
        y.objective_coefficient = -0.5
        y.lower_bound = 0.
        x.add_metabolites({c: 1})
        y.add_metabolites({c: 1})
        m = Model()
        m.add_reactions([x, y])
        lp = self.solver.create_problem(m)
        quadratic_obj = scipy.sparse.eye(2) * 2
        solver.set_quadratic_objective(lp, quadratic_obj)
        solver.solve_problem(lp, objective_sense="minimize")
        solution = solver.format_solution(lp, m)
        self.assertEqual(solution.status, "optimal")
        # Respecting linear objectives also makes the objective value 1.
        self.assertAlmostEqual(solution.f, 1.)
        self.assertAlmostEqual(solution.x_dict["y"], 1.)
        self.assertAlmostEqual(solution.x_dict["y"], 1.)
        # When the linear objectives are removed the objective value is 2.
        solver.change_variable_objective(lp, 0, 0.)
        solver.change_variable_objective(lp, 1, 0.)
        solver.solve_problem(lp, objective_sense="minimize")
        solution = solver.format_solution(lp, m)
        self.assertEqual(solution.status, "optimal")
        self.assertAlmostEqual(solution.f, 2.)
        # test quadratic from solve function
        solution = solver.solve(m, quadratic_component=quadratic_obj,
                                objective_sense="minimize")
        self.assertEqual(solution.status, "optimal")
        self.assertAlmostEqual(solution.f, 1.)
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
        self.assertEqual(solution.status, "optimal")
        self.assertAlmostEqual(solution.f, 6)
        self.assertAlmostEqual(solution.x_dict["x"], 2, places=6)
        self.assertAlmostEqual(solution.x_dict["y"], 2, places=6)
        self.assertAlmostEqual(solution.x_dict["z"], 2, places=6)

for solver_name in solvers.solver_dict:
    exec('class %sTester(TestCobraSolver, TestCase): None' % solver_name)
    exec('%sTester.solver_name = "%s"' % (solver_name, solver_name))

# make a test suite to run all of the tests
loader = TestLoader()
suite = loader.loadTestsFromModule(sys.modules[__name__])


def test_all():
    TextTestRunner(verbosity=2).run(suite)

if __name__ == "__main__":
    test_all()
