import unittest
import warnings
import os
try:
    from cPickle import load
except:
    from pickle import load
import sys

from cobra.manipulation import initialize_growth_medium
from cobra.test import create_test_model
from cobra import Model, Reaction, Metabolite
from cobra import solvers
from cobra.solvers import __legacy_solver

class TestCobraSolver(unittest.TestCase):
    def setUp(self):
        self.model = create_test_model()
        initialize_growth_medium(self.model, 'MgM')
        self.old_solution = 0.320064
        self.infeasible_problem = Model()
        metabolite_1 = Metabolite("met1")
        metabolite_2 = Metabolite("met2")
        reaction_1 = Reaction("rxn1")
        reaction_2 = Reaction("rxn2")
        reaction_1.add_metabolites({metabolite_1: 1})
        reaction_2.add_metabolites({metabolite_2: 1})
        reaction_1.lower_bound = 1
        reaction_2.upper_bound = 2
        self.infeasible_problem.add_reactions([reaction_1, reaction_2])
        #self.infeasible_problem.update()


def add_new_test(TestCobraSolver, solver_name, solver):
    """Creates a test set for each of the solvers that are installed
    using the modular interface.

    
    """
    def test_attributes(self):
        self.assertTrue(hasattr(solver, "create_problem"))
        self.assertTrue(hasattr(solver, "solve_problem"))
        # self.assertTrue(hasattr(solver, "update_problem"))
    def test_setup(self):
        solver.create_problem(self.model)
    def test_solve_feasible(self):
        solution = solver.solve(self.model)
        self.assertEqual(solution.status, "optimal")
        self.assertAlmostEqual(self.old_solution, \
            solution.objective_value, places=4)
    def test_solve_infeasible(self):
        solution = solver.solve(self.infeasible_problem)
        self.assertEqual(solution.status, "infeasible")
    def test_independent_creation(self):
        feasible_lp = solver.create_problem(self.model)
        infeasible_lp = solver.create_problem(self.infeasible_problem)
        feasible_solution = solve_problem(lp)
        infeasible_solution = solve_problem(lp)
        self.assertEqual(feasible_solution.status, "optimal")
        self.assertAlmostEqual(self.old_solution, \
            feasible_solution.objective_value, places=4)
        self.assertEqual(infeasible_solution.status, "infeasible")
    setattr(TestCobraSolver, "test_%s_create" % solver_name, \
        test_setup)
    setattr(TestCobraSolver, "test_%s_attributes" % solver_name, \
        test_attributes)
    setattr(TestCobraSolver, "test_%s_feasible_solve" % solver_name, \
        test_solve_feasible)
    setattr(TestCobraSolver, "test_%s_infeasible_solve" % solver_name, \
        test_solve_infeasible)
    setattr(TestCobraSolver, "test_%s_independent_creation" % solver_name, \
        test_solve_infeasible)

def add_legacy_test(TestCobraSolver, solver_name, solver_function):
    """Creates a test set for each of the installed solvers using the
    legacy interface.

    """
    def test_solve_feasible(self):

        the_solution = solver_function(self.model)['the_solution']
        self.assertEqual(the_solution.status, 'optimal')
        self.assertAlmostEqual(self.old_solution, the_solution.f, places=4)
    setattr(TestCobraSolver, "test_%s_feasible_solve" % solver_name, \
        test_solve_feasible)
        


if not __legacy_solver:
    add_test = add_new_test
else:
    add_test = add_legacy_test

    
for solver_name, solver in solvers.solver_dict.iteritems():
    add_test(TestCobraSolver, solver_name, solver)
# make a test suite to run all of the tests
loader = unittest.TestLoader()
suite = loader.loadTestsFromModule(sys.modules[__name__])

def test_all():
    unittest.TextTestRunner(verbosity=2).run(suite)

if __name__ == "__main__":
    test_all()
