from unittest import TestCase, TestLoader, TextTestRunner, skipIf
from warnings import catch_warnings
import os, sys
sys.path.insert(0, "../..")
from cobra.test import data_directory, create_test_model
from cobra.test import salmonella_sbml as test_sbml_file
from cobra.test import salmonella_pickle as test_pickle
from cobra import Object, Model, Metabolite, Reaction, io, DictList
sys.path.pop(0)
#from .. import flux_analysis

# libraries which may or may not be installed
libraries = ["libsbml", "glpk", "gurobipy", "cplex"]
for library in libraries:
    try:
        exec("import %s" % library)
    except ImportError:
        exec("%s = None" % library)


class TestDictList(TestCase):
    def setUp(self):
        self.obj = Object("test1")
        self.list = DictList()
        self.list.append(self.obj)

    def testAppend(self):
        obj2 = Object("test2")
        self.list.append(obj2)
        self.assertRaises(ValueError, self.list.append,
            Object("test1"))
        self.assertEqual(self.list.index(obj2), 1)
        self.assertEqual(self.list[1], obj2)
        self.assertEqual(self.list.get_by_id("test2"), obj2)
        self.assertEqual(len(self.list), 2)

    def testExtend(self):
        obj_list = [Object("test%d" % (i)) for i in range(2, 10)]
        self.list.extend(obj_list)
        self.assertEqual(self.list[1].id, "test2")
        self.assertEqual(self.list.get_by_id("test2"), obj_list[0])
        self.assertEqual(self.list[8].id, "test9")
        self.assertEqual(len(self.list), 9)

    def testIadd(self):
        obj_list = [Object("test%d" % (i)) for i in range(2, 10)]
        self.list += obj_list
        self.assertEqual(self.list[1].id, "test2")
        self.assertEqual(self.list.get_by_id("test2"), obj_list[0])
        self.assertEqual(self.list[8].id, "test9")
        self.assertEqual(len(self.list), 9)

    def testAdd(self):
        obj_list = [Object("test%d" % (i)) for i in range(2, 10)]
        sum = self.list + obj_list
        self.assertEqual(self.list[0].id, "test1")
        self.assertEqual(sum[1].id, "test2")
        self.assertEqual(sum.get_by_id("test2"), obj_list[0])
        self.assertEqual(sum[8].id, "test9")
        self.assertEqual(len(self.list), 1)
        self.assertEqual(len(sum), 9)


class CobraTestCase(TestCase):
    def setUp(self):
        self.model = create_test_model(test_pickle)


class TestCobraCore(CobraTestCase):
    """test core cobra functions"""

    def test_add_reaction(self):
        """test adding and deleting reactions"""
        old_reaction_count = len(self.model.reactions)
        old_metabolite_count = len(self.model.metabolites)
        dummy_metabolite_1 = Metabolite("test_foo_1")
        dummy_metabolite_2 = Metabolite("test_foo_2")
        dummy_reaction = Reaction("test_foo_reaction")
        dummy_reaction.add_metabolites({dummy_metabolite_1: -1,
                                        dummy_metabolite_2: 1})
        self.model.add_reaction(dummy_reaction)
        self.assertEqual(self.model.reactions._object_dict[dummy_reaction.id],
                         dummy_reaction)
        for x in [dummy_metabolite_1, dummy_metabolite_2]:
            self.assertEqual(self.model.metabolites._object_dict[x.id],
                             x)
            

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


class TestCobraIO(CobraTestCase):

    @skipIf(libsbml is None, "libsbml is required")
    def test_sbml_read(self):
        with catch_warnings(record=True) as w:
            model = io.read_sbml_model(test_sbml_file)
        self.assertEqual(len(model.reactions), len(self.model.reactions))
        # make sure that an error is raised when given a nonexistent file
        self.assertRaises(IOError, io.read_sbml_model,
            "fake_file_which_does_not_exist")

    @skipIf(libsbml is None, "libsbml is required")
    def test_sbml_write(self):
        io.write_sbml_model(self.model, "test_sbml_write.xml")

# make a test suite to run all of the tests
loader = TestLoader()
suite = loader.loadTestsFromModule(sys.modules[__name__])

def test_all():
    TextTestRunner(verbosity=2).run(suite)

if __name__ == "__main__":
    test_all()
