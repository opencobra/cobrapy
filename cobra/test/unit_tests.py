from __future__ import with_statement
import sys
from warnings import warn  # TODO - catch known warnings
from unittest import TestCase, TestLoader, TextTestRunner
from tempfile import gettempdir
from os import unlink
from copy import deepcopy
from os.path import join
try:  #skipIf is not in python 2.6 / 2.5, so use unittest2
    from unittest import skipIf
except:
    from unittest2 import skipIf


# deal with absolute imports by adding the appropriate directory to the path
if __name__ == "__main__":
    sys.path.insert(0, "../..")
    from cobra.test import data_directory, create_test_model
    from cobra.test import salmonella_sbml as test_sbml_file
    from cobra.test import salmonella_pickle as test_pickle
    from cobra import Object, Model, Metabolite, Reaction, io, DictList
    sys.path.pop(0)
    #assert 0
else:
    from . import data_directory, create_test_model
    from . import salmonella_sbml as test_sbml_file
    from . import salmonella_pickle as test_pickle
    from .. import Object, Model, Metabolite, Reaction, io, DictList

# libraries which may or may not be installed
libraries = ["glpk", "gurobipy", "cplex"]
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

    def testDeepcopy(self):
        from copy import deepcopy
        copied = deepcopy(self.list)
        for i, v in enumerate(self.list):
            assert self.list[i].id == copied[i].id
            assert self.list[i] is not copied[i]

    def testQuery(self):
        obj2 = Object("test2")
        self.list.append(obj2)
        result = self.list.query("test1")  # matches only test1
        self.assertEqual(len(result), 1)
        self.assertEqual(result[0], self.obj)
        result = self.list.query("test")  # matches test1 and test2
        self.assertEqual(len(result), 2)

            

class CobraTestCase(TestCase):
    def setUp(self):
        self.model = create_test_model(test_pickle)


class TestCobraCore(CobraTestCase):
    """test core cobra functions"""

    def test_add_reaction(self):
        old_reaction_count = len(self.model.reactions)
        old_metabolite_count = len(self.model.metabolites)
        dummy_metabolite_1 = Metabolite("test_foo_1")
        dummy_metabolite_2 = Metabolite("test_foo_2")
        actual_metabolite = self.model.metabolites[0]
        copy_metabolite = self.model.metabolites[1].copy()
        dummy_reaction = Reaction("test_foo_reaction")
        dummy_reaction.add_metabolites({dummy_metabolite_1: -1,
                                        dummy_metabolite_2: 1,
                                        copy_metabolite: -2,
                                        actual_metabolite: 1})
        self.model.add_reaction(dummy_reaction)
        self.assertEqual(self.model.reactions.get_by_id(dummy_reaction.id),
                         dummy_reaction)
        for x in [dummy_metabolite_1, dummy_metabolite_2]:
            self.assertEqual(self.model.metabolites.get_by_id(x.id), x)
        # should have added 1 reaction and 2 metabolites
        self.assertEqual(len(self.model.reactions), old_reaction_count + 1)
        self.assertEqual(len(self.model.metabolites), old_metabolite_count + 2)
        # tests on theadded reaction
        reaction_in_model = self.model.reactions.get_by_id(dummy_reaction.id)
        self.assertIs(type(reaction_in_model), Reaction)
        self.assertIs(reaction_in_model, dummy_reaction)
        self.assertEqual(len(reaction_in_model._metabolites), 4)
        for i in reaction_in_model._metabolites:
            self.assertEqual(type(i), Metabolite)
        # tests on the added metabolites
        met1_in_model = self.model.metabolites.get_by_id(dummy_metabolite_1.id)
        self.assertIs(met1_in_model, dummy_metabolite_1)
        #assertIsNot is not in python 2.6
        copy_in_model = self.model.metabolites.get_by_id(copy_metabolite.id)
        self.assertTrue(copy_metabolite is not copy_in_model)
        self.assertIs(type(copy_in_model), Metabolite)
        self.assertTrue(dummy_reaction in actual_metabolite._reaction)

    def test_delete_reaction(self):
        old_reaction_count = len(self.model.reactions)
        self.model.remove_reactions([self.model.reactions.get_by_id("PGI")])
        self.assertEqual(len(self.model.reactions), old_reaction_count - 1)
        with self.assertRaises(KeyError):
            self.model.reactions.get_by_id("PGI")
        # TODO - delete by id - will this be supported?
        # TODO - delete orphan metabolites - will this be expected behavior?
        

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


    def test_deepcopy(self):
        """Verify that reference structures are maintained when deepcopying.
        
        """
        model_copy = deepcopy(self.model)
        for gene, gene_copy in zip(self.model.genes, model_copy.genes):
            self.assertEqual(gene.id, gene_copy.id)
            reactions = list(gene.get_reaction())
            reactions.sort()
            reactions_copy = list(gene_copy.get_reaction())
            reactions_copy.sort()
            self.assertEqual(reactions, reactions_copy)
        for reaction, reaction_copy in zip(self.model.reactions, model_copy.reactions):
            self.assertEqual(reaction.id, reaction_copy.id)
            metabolites = reaction._metabolites.keys()
            metabolites.sort()
            metabolites_copy = reaction_copy._metabolites.keys()
            metabolites_copy.sort()
            self.assertEqual(metabolites, metabolites_copy)
        
        

class TestCobraIO(CobraTestCase):
    try:
        from cobra.io import sbml
        __test_sbml = True
    except:
        __test_sbml = False
    try:
        from scipy.io import loadmat
        __test_matlab = True
    except:
        __test_matlab = False

    @skipIf(not __test_sbml, "libsbml required")
    def test_sbml_read(self):
        from warnings import catch_warnings
        with catch_warnings(record=True) as w:
            model = io.read_sbml_model(test_sbml_file)
        self.assertEqual(len(model.reactions), len(self.model.reactions))
        # make sure that an error is raised when given a nonexistent file
        self.assertRaises(IOError, io.read_sbml_model,
                          "fake_file_which_does_not_exist")

    @skipIf(not __test_sbml, "libsbml required")
    def test_sbml_write(self):
        test_output_filename = join(gettempdir(), 'test_sbml_write.xml')
        io.write_sbml_model(self.model, test_output_filename)
        #cleanup the test file
        unlink(test_output_filename)
    
    @skipIf(not __test_matlab, "scipy.io.loadmat required")
    def test_mat_read_write(self):
        from warnings import catch_warnings
        test_output_filename = join(gettempdir(), "test_mat_write.mat")
        io.save_matlab_model(self.model, test_output_filename)
        with catch_warnings(record=True) as w:
            reread = io.load_matlab_model(test_output_filename)
        self.assertEqual(len(self.model.reactions), len(reread.reactions))
        self.assertEqual(len(self.model.metabolites), len(reread.metabolites))
        for i in range(len(self.model.reactions)):
            self.assertEqual(len(self.model.reactions[i]._metabolites), \
                len(reread.reactions[i]._metabolites))
            self.assertEqual(self.model.reactions[i].id, reread.reactions[i].id)
        unlink(test_output_filename)

# make a test suite to run all of the tests
loader = TestLoader()
suite = loader.loadTestsFromModule(sys.modules[__name__])

def test_all():
    TextTestRunner(verbosity=2).run(suite)

if __name__ == "__main__":
    test_all()
