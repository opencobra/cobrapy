import unittest
import warnings
import os
try:
    from cPickle import load
except:
    from pickle import load
import sys


# from . import data_directory, ecoli_sbml, ecoli_pickle, create_test_model
from cobra.test import create_test_model
from cobra.test import salmonella_sbml as test_sbml_file
from cobra.test import salmonella_pickle as test_pickle
from cobra import Model, Reaction, Metabolite
from cobra import solvers


class TestCobraSolver(unittest.TestCase):
    def setUp(self):
        self.model = create_test_model()
        self.old_solution = 0.982371812727
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
        self.infeasible_problem.update()


def add_test(TestCobraSolver, solver_name, solver):
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

for solver_name, solver in solvers.solver_dict.iteritems():
    add_test(TestCobraSolver, solver_name, solver)
# make a test suite to run all of the tests
loader = unittest.TestLoader()
suite = loader.loadTestsFromModule(sys.modules[__name__])

def test_all():
    unittest.TextTestRunner(verbosity=2).run(suite)

if __name__ == "__main__":
    test_all()
