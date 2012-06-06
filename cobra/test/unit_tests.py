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


class TestDictList(unittest.TestCase):
    def setUp(self):
        self.obj = cobra.Object.Object("test1")
        self.list = cobra.collections.DictList.DictList()
        self.list.append(self.obj)

    def testAppend(self):
        obj2 = cobra.Object.Object("test2")
        self.list.append(obj2)
        self.assertRaises(ValueError, self.list.append,
            cobra.Object.Object("test1"))
        self.assertEqual(self.list.index(obj2), 1)
        self.assertEqual(self.list[1], obj2)
        self.assertEqual(self.list.get_by_id("test2"), obj2)
        self.assertEqual(len(self.list), 2)

    def testExtend(self):
        obj_list = [cobra.Object.Object("test%d" % (i)) for i in range(2, 10)]
        self.list.extend(obj_list)
        self.assertEqual(self.list[1].id, "test2")
        self.assertEqual(self.list.get_by_id("test2"), obj_list[0])
        self.assertEqual(self.list[8].id, "test9")
        self.assertEqual(len(self.list), 9)

    def testIadd(self):
        obj_list = [cobra.Object.Object("test%d" % (i)) for i in range(2, 10)]
        self.list += obj_list
        self.assertEqual(self.list[1].id, "test2")
        self.assertEqual(self.list.get_by_id("test2"), obj_list[0])
        self.assertEqual(self.list[8].id, "test9")
        self.assertEqual(len(self.list), 9)

    def testAdd(self):
        obj_list = [cobra.Object.Object("test%d" % (i)) for i in range(2, 10)]
        sum = self.list + obj_list
        self.assertEqual(self.list[0].id, "test1")
        self.assertEqual(sum[1].id, "test2")
        self.assertEqual(sum.get_by_id("test2"), obj_list[0])
        self.assertEqual(sum[8].id, "test9")
        self.assertEqual(len(self.list), 1)
        self.assertEqual(len(sum), 9)


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


def test_optimizer(FluxTest, solver):
    """function to perform unittests on a generic solver"""
    old_solution = FluxTest.model.solution.f
    opt_function = cobra.flux_analysis.solvers.__getattribute__(
        "optimize_%s" % (solver))
    # test a feasible model
    # using model.optimize
    FluxTest.model.optimize(solver=solver)
    FluxTest.assertEqual(FluxTest.model.solution.status, "optimal")
    FluxTest.assertAlmostEqual(FluxTest.model.solution.f,
        old_solution, places=5)
    # using the optimization function directly
    solution = opt_function(FluxTest.model)["the_solution"]
    FluxTest.assertEqual(solution.status, "optimal")
    FluxTest.assertAlmostEqual(solution.f, old_solution, places=5)
    # test an infeasible model
    # using model.optimize
    FluxTest.infeasible_problem.optimize(solver=solver)
    FluxTest.assertEqual(FluxTest.infeasible_problem.solution.status,
        "infeasible")
    # using the optimization function directly
    solution = opt_function(FluxTest.infeasible_problem)["the_solution"]
    FluxTest.assertEqual(solution.status, "infeasible")

    # ensure a copied model can also be solved
    model_copy = FluxTest.model.copy()
    model_copy.optimize(solver=solver)
    FluxTest.assertAlmostEqual(model_copy.solution.f, old_solution, places=5)
    solution = opt_function(model_copy)["the_solution"]
    FluxTest.assertAlmostEqual(solution.f, old_solution, places=5)


class TestCobraFlux_analysis(CobraTestCase):

    def setUp(self):
        CobraTestCase.setUp(self)
        self.old_solution = self.model.solution.f
        self.infeasible_problem = cobra.Model()
        metabolite_1 = cobra.Metabolite("met1")
        metabolite_2 = cobra.Metabolite("met2")
        reaction_1 = cobra.Reaction("rxn1")
        reaction_2 = cobra.Reaction("rxn2")
        reaction_1.add_metabolites({metabolite_1: 1})
        reaction_2.add_metabolites({metabolite_2: 1})
        reaction_1.lower_bound = 1
        reaction_2.upper_bound = 2
        self.infeasible_problem.add_reactions([reaction_1, reaction_2])

    @unittest.skipIf(glpk is None, "glpk is required")
    def test_glpk_optimization(self):
        test_optimizer(self, "glpk")

    @unittest.skipIf(gurobipy is None, "gurobipy is required")
    def test_gurobi_optimization(self):
        test_optimizer(self, "gurobi")

    @unittest.skipIf(cplex is None, "cplex is required")
    def test_cplex_optimization(self):
        test_optimizer(self, "cplex")


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
