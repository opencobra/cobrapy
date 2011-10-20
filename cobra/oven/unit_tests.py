import unittest
import warnings
from cPickle import load
import sys

import cobra
from cobra.test import data_directory, salmonella_sbml, salmonella_pickle

# libraries which may or may not be installed
libraries = ["libsbml", "glpk", "gurobipy", "cplex"]
for library in libraries:
    try:
        exec("import %s" % library)
    except ImportError:
        exec("%s = None" % library)


class CobraTestCase(unittest.TestCase):
    def setUp(self):
        infile = open(salmonella_pickle, "rb")
        self.model = load(infile)
        infile.close()


class TestCobraCore(CobraTestCase):
    """test core cobra functions"""

    def test_add_reaction(self):
        """test adding and deleting reactions"""
        old_reaction_count = len(self.model.reactions)
        old_metabolite_count = len(self.model.metabolites)
        dummy_metabolite_1 = cobra.Metabolite("test_foo_1")
        dummy_metabolite_2 = cobra.Metabolite("test_foo_2")
        dummy_reaction = cobra.Reaction("test_foo_reaction")
        dummy_reaction.add_metabolites({dummy_metabolite_1: -1,
                                        dummy_metabolite_2: 1})
        self.model.add_reaction(dummy_reaction)
        self.model.update()
        self.assertEqual(len(self.model.reactions), self.model._S.shape[1])
        self.assertEqual(len(self.model.metabolites), self.model._S.shape[0])
        self.assertEqual(len(self.model.reactions), old_reaction_count + 1)
        self.assertEqual(len(self.model.metabolites),
            old_metabolite_count + 2)
        self.model.remove_reactions(dummy_reaction)
        self.assertEqual(len(self.model.reactions), old_reaction_count)
        self.model.update()
        self.assertEqual(len(self.model.reactions), self.model._S.shape[1])
        self.assertEqual(len(self.model.metabolites), self.model._S.shape[0])

    def test_copy(self):
        """modifying copy should not modify the original"""
        # test that deleting reactions in the copy does not change the
        # number of reactions in the original model
        model_copy = self.model.copy()
        old_reaction_count = len(self.model.reactions)
        self.assertEqual(len(self.model.reactions), len(model_copy.reactions))
        self.assertEqual(len(self.model.metabolites),
            len(model_copy.metabolites))
        model_copy.remove_reactions(model_copy.reactions[0:5])
        self.assertEqual(old_reaction_count, len(self.model.reactions))
        self.assertNotEqual(len(self.model.reactions),
            len(model_copy.reactions))


class TestCobraFlux_analysis(CobraTestCase):

    def setUp(self):
        CobraTestCase.setUp(self)
        self.old_solution = self.model.solution.f

    @unittest.skipIf(glpk is None, "glpk is required")
    def test_glpk_optimization(self):
        old_solution = self.model.solution.f
        self.model.optimize(solver="glpk")
        self.assertAlmostEqual(self.model.solution.f, old_solution, places=5)

    @unittest.skipIf(gurobipy is None, "gurobipy is required")
    def test_gurobi_optimization(self):
        old_solution = self.model.solution.f
        self.model.optimize(solver="gurobi")
        self.assertAlmostEqual(self.model.solution.f, old_solution, places=5)

    @unittest.skipIf(cplex is None, "cplex is required")
    def test_cplex_optimization(self):
        old_solution = self.model.solution.f
        self.model.optimize(solver="cplex")
        self.assertAlmostEqual(self.model.solution.f, old_solution, places=5)


class TestCobraIO(CobraTestCase):

    @unittest.skipIf(libsbml is None, "libsbml is required")
    def test_sbml_read(self):
        with warnings.catch_warnings(record=True) as w:
            from cobra.io.sbml import create_cobra_model_from_sbml_file
            model = create_cobra_model_from_sbml_file(salmonella_sbml)
        self.assertEqual(len(model.reactions), len(self.model.reactions))
        # make sure that an error is raised when given a nonexistent file
        self.assertRaises(IOError, create_cobra_model_from_sbml_file,
            "not_a_real_file_at_all")

    @unittest.skipIf(libsbml is None, "libsbml is required")
    def test_sbml_write(self):
        from cobra.io.sbml import write_cobra_model_to_sbml_file
        write_cobra_model_to_sbml_file(self.model, "salmonella_out.xml")

# make a test suite to run all of the tests
loader = unittest.TestLoader()
suite = loader.loadTestsFromModule(sys.modules[__name__])

if __name__ == "__main__":
    unittest.TextTestRunner(verbosity=2).run(suite)
