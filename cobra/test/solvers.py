from unittest import TestCase, TestLoader, TextTestRunner
import sys
# deal with absolute imports by adding the appropriate directory to the path
if __name__ == "__main__":
    sys.path.insert(0, "../..")
    from cobra.manipulation import initialize_growth_medium
    from cobra.test import create_test_model
    from cobra import Model, Reaction, Metabolite
    from cobra import solvers
    sys.path.pop(0)  # remove the added directory to the path
else:
    from ..manipulation import initialize_growth_medium
    from . import create_test_model
    from .. import Model, Reaction, Metabolite
    from .. import solvers

solver_dict = solvers.solver_dict

class TestCobraSolver(TestCase):
    def setUp(self):
        self.model = create_test_model()
        initialize_growth_medium(self.model, 'MgM')
        self.old_solution = 0.320064
        self.infeasible_model = Model()
        metabolite_1 = Metabolite("met1")
        #metabolite_2 = Metabolite("met2")
        reaction_1 = Reaction("rxn1")
        reaction_2 = Reaction("rxn2")
        reaction_1.add_metabolites({metabolite_1: 1})
        reaction_2.add_metabolites({metabolite_1: 1})
        reaction_1.lower_bound = 1
        reaction_2.upper_bound = 2
        self.infeasible_model.add_reactions([reaction_1, reaction_2])
        #self.infeasible_model.update()


def add_new_test(TestCobraSolver, solver_name, solver):
    """Creates a test set for each of the solvers that are installed
    using the modular interface.

    """
    def attributes(self):
        solver = solver_dict[solver_name]
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

    def creation(self):
        solver = solver_dict[solver_name]
        solver.create_problem(self.model)

    def solve_feasible(self):
        solver = solver_dict[solver_name]
        solution = solver.solve(self.model)
        self.assertEqual(solution.status, "optimal")
        self.assertAlmostEqual(self.old_solution, \
            solution.f, places=4)

    def solve_minimize(self):
        solver = solver_dict[solver_name]
        solution = solver.solve(self.model, objective_sense='minimize')
        self.assertEqual(solution.status, "optimal")
        self.assertAlmostEqual(0, solution.f, places=4)

    def low_level_control(self):
        solver = solver_dict[solver_name]
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


    def set_objective_sense(self):
        maximize = solver.create_problem(self.model, objective_sense="maximize")
        minimize = solver.create_problem(self.model, objective_sense="minimize")
        solver.solve_problem(maximize)
        solver.solve_problem(minimize)
        max_solution = solver.format_solution(maximize, self.model)
        min_solution = solver.format_solution(minimize, self.model)
        self.assertAlmostEqual(0, min_solution.f, places=4)
        self.assertAlmostEqual(self.old_solution, max_solution.f, places=4)
        # if we set minimize at creation, can we override it at solve
        solver.solve_problem(minimize, objective_sense="maximize")
        override_minimize = solver.format_solution(minimize, self.model)
        self.assertAlmostEqual(max_solution.f, override_minimize.f, places=4)

    def solve_mip(self):
        solver = solver_dict[solver_name]
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
        self.assertAlmostEqual(int_sol.f, 2.2)
        self.assertAlmostEqual(int_sol.x_dict["y"], 2.0)

    def solve_infeasible(self):
        solution = solver.solve(self.infeasible_model)
        self.assertEqual(solution.status, "infeasible")

    def independent_creation(self):
        feasible_lp = solver.create_problem(self.model)
        infeasible_lp = solver.create_problem(self.infeasible_model)
        solver.solve_problem(feasible_lp)
        solver.solve_problem(infeasible_lp)
        feasible_solution = solver.format_solution(feasible_lp, self.model)
        infeasible_solution = solver.format_solution(infeasible_lp, self.infeasible_model)
        self.assertEqual(feasible_solution.status, "optimal")
        self.assertAlmostEqual(self.old_solution, \
            feasible_solution.f, places=4)
        self.assertEqual(infeasible_solution.status, "infeasible")

    for tester in [attributes, creation, solve_feasible, solve_minimize,
                    set_objective_sense, solve_mip, solve_infeasible,
                    independent_creation, low_level_control]:
        test_name = tester.__name__ if hasattr(tester, "__name__") \
            else tester.func_name
        setattr(TestCobraSolver, "test_%s_%s" %
                (solver_name, test_name), tester)


for solver_name, solver in solvers.solver_dict.items():
    add_new_test(TestCobraSolver, solver_name, solver)
# make a test suite to run all of the tests
loader = TestLoader()
suite = loader.loadTestsFromModule(sys.modules[__name__])

def test_all():
    TextTestRunner(verbosity=2).run(suite)

if __name__ == "__main__":
    test_all()
