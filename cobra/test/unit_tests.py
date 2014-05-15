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
    from cobra.flux_analysis import single_deletion
    sys.path.pop(0)
    #assert 0
else:
    from . import data_directory, create_test_model
    from . import salmonella_sbml as test_sbml_file
    from . import salmonella_pickle as test_pickle
    from .. import Object, Model, Metabolite, Reaction, io, DictList
    from ..flux_analysis import single_deletion
# libraries which may or may not be installed
libraries = ["glpk", "gurobipy", "cplex", "scipy"]
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
            self.assertEqual(self.list[i].id, copied[i].id)
            self.assertIsNot(self.list[i], copied[i])

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
        self.model_class = Model


class TestReactions(CobraTestCase):
    def testGPR(self):
        model = self.model_class()
        reaction = Reaction("test")
        # set a gpr to  reaction not in a model
        reaction.gene_reaction_rule = "(g1 or g2) and g3"
        self.assertEqual(reaction.gene_reaction_rule, "(g1 or g2) and g3")
        self.assertEqual(len(reaction.genes), 3)
        # adding reaction with a GPR propagates to the model
        model.add_reaction(reaction)
        self.assertEqual(len(model.genes), 3)
        # ensure the gene objects are the same in the model and reaction
        reaction_gene = list(reaction.genes)[0]
        model_gene = model.genes.get_by_id(reaction_gene.id)
        self.assertIs(reaction_gene, model_gene)


    def testGPR_modification(self):
        model = self.model
        reaction = model.reactions.get_by_id("PGI")
        old_gene = list(reaction.genes)[0]
        old_gene_reaction_rule = reaction.gene_reaction_rule
        new_gene = model.genes.get_by_id("s0001")
        # add an existing 'gene' to the gpr
        reaction.gene_reaction_rule = 's0001'
        self.assertIn(new_gene, reaction.genes)
        self.assertIn(reaction, new_gene.reactions)
        # removed old gene correctly
        self.assertNotIn(old_gene, reaction.genes)
        self.assertNotIn(reaction, old_gene.reactions)
        # add a new 'gene' to the gpr
        reaction.gene_reaction_rule = 'fake_gene'
        self.assertTrue(model.genes.has_id("fake_gene"))
        fake_gene = model.genes.get_by_id("fake_gene")
        self.assertIn(fake_gene, reaction.genes)
        self.assertIn(reaction, fake_gene.reactions)

    def test_add_metabolite(self):
        """adding a metabolite to a reaction in a model"""
        model = self.model
        reaction = model.reactions.get_by_id("PGI")
        reaction.add_metabolites({model.metabolites[0]: 1})
        self.assertIn(model.metabolites[0], reaction.metabolites)
        fake_metabolite = Metabolite("fake")
        reaction.add_metabolites({fake_metabolite: 1})
        self.assertIn(fake_metabolite, reaction.metabolites)
        self.assertTrue(model.metabolites.has_id("fake"))
        self.assertIs(model.metabolites.get_by_id("fake"), fake_metabolite)

class TestCobraModel(CobraTestCase):
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

    def test_add_reaction_from_other_model(self):
        model = self.model
        other = model.copy()
        for i in other.reactions:
            i.id += "_other"
        other.reactions._generate_index()
        model.add_reactions(other.reactions)

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
        """Reference structures are maintained when deepcopying"""
        model_copy = deepcopy(self.model)
        for gene, gene_copy in zip(self.model.genes, model_copy.genes):
            self.assertEqual(gene.id, gene_copy.id)
            reactions = sorted(i.id for i in gene.reactions)
            reactions_copy = sorted(i.id for i in gene_copy.reactions)
            self.assertEqual(reactions, reactions_copy)
        for reaction, reaction_copy in zip(self.model.reactions, model_copy.reactions):
            self.assertEqual(reaction.id, reaction_copy.id)
            metabolites = sorted(i.id for i in reaction._metabolites)
            metabolites_copy = sorted(i.id for i in reaction_copy._metabolites)
            self.assertEqual(metabolites, metabolites_copy)

    def test_add_reaction(self):
        """test reaction addition

        Need to verify that no orphan genes or metabolites are
        contained in reactions after adding them to the model.
        """
        _model = self.model_class('test')
        _model.add_reactions([x.copy() for x in self.model.reactions])
        _genes = []
        _metabolites = []
        [(_genes.extend(x.genes), _metabolites.extend(x.metabolites))
         for x in _model.reactions];
        _orphan_genes = [x for x in _genes if x.model is not _model]
        _orphan_metabolites = [x for x in _metabolites if x.model is not _model]
        self.assertEqual(len(_orphan_genes), 0, msg='It looks like there are dangling genes when running Model.add_reactions')
        self.assertEqual(len(_orphan_metabolites), 0, msg='It looks like there are dangling metabolites when running Model.add_reactions')

@skipIf(scipy is None, "scipy required for ArrayBasedModel")
class TestCobraArrayModel(TestCobraModel):
    def setUp(self):
        model = create_test_model(test_pickle).to_array_based_model()
        self.model_class = model.__class__
        self.model = model

    def test_array_based_model(self):
        for matrix_type in ["scipy.dok_matrix", "scipy.lil_matrix"]:
            model = create_test_model().to_array_based_model(matrix_type=matrix_type)
            self.assertEqual(model.S[0, 0], -1)
            self.assertEqual(model.S[43, 0], 0)
            model.S[43, 0] = 1
            self.assertEqual(model.S[43, 0], 1)
            self.assertEqual(model.reactions[0].metabolites[model.metabolites[43]], 1)
            model.S[43, 0] = 0
            self.assertEqual(model.lower_bounds[0], model.reactions[0].lower_bound)
            self.assertEqual(model.lower_bounds[5], model.reactions[5].lower_bound)
            self.assertEqual(model.upper_bounds[0], model.reactions[0].upper_bound)
            self.assertEqual(model.upper_bounds[5], model.reactions[5].upper_bound)
            model.lower_bounds[6] = 2
            self.assertEqual(model.lower_bounds[6], 2)
            self.assertEqual(model.reactions[6].lower_bound, 2)
            # this should fail because it is the wrong size
            with self.assertRaises(Exception):
                model.upper_bounds = [0, 1]
            model.upper_bounds = [0] * len(model.reactions)
            self.assertEqual(max(model.upper_bounds), 0)

    def test_array_based_model_add(self):
        for matrix_type in ["scipy.dok_matrix", "scipy.lil_matrix"]:
            model = create_test_model().to_array_based_model(matrix_type=matrix_type)
            test_reaction = Reaction("test")
            test_reaction.add_metabolites({model.metabolites[0]: 4})
            test_reaction.lower_bound = -3.14
            model.add_reaction(test_reaction)
            self.assertEqual(len(model.reactions), 2547)
            self.assertEqual(model.S.shape[1], 2547)
            self.assertEqual(len(model.lower_bounds), 2547)
            self.assertEqual(model.S[0, 2546], 4)
            self.assertEqual(model.S[0, 0], -1)
            self.assertEqual(model.lower_bounds[2546], -3.14)

class TestCobraIO(CobraTestCase):
    try:
        from cobra.io import sbml
        __test_sbml = True
    except:
        __test_sbml = False
    try:
        from libsbml import FbcExtension
        __test_sbml_fbc = True
    except:
        __test_sbml_fbc = False
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
        self.assertEqual(len(model.metabolites), len(self.model.metabolites))
        # make sure that an error is raised when given a nonexistent file
        self.assertRaises(IOError, io.read_sbml_model,
                          "fake_file_which_does_not_exist")

    @skipIf(not __test_sbml, "libsbml required")
    def test_sbml_write(self):
        test_output_filename = join(gettempdir(), 'test_sbml_write.xml')
        io.write_sbml_model(self.model, test_output_filename, use_fbc_package=False)
        #cleanup the test file
        unlink(test_output_filename)

    @skipIf(not __test_sbml_fbc, "libsbml with fbc package required")
    def test_sbml_fbc(self):
        """This tests whether activating the fbc extensions affect simulation results.
        
        """
        test_output_filename = join(gettempdir(), 'test_sbml_write.xml')
        #Save the model in SBML+fbc
        io.write_sbml_model(self.model, test_output_filename)
        fbc_model = io.read_sbml_model(test_output_filename)
        #cleanup the test file
        unlink(test_output_filename)
        #Compare initial optimizations
        fbc_model.optimize()
        self.model.optimize()
        self.assertAlmostEqual(fbc_model.solution.f, self.model.solution.f, places=3)
        
        #Compare the single deletion results of the first 100 genes
        gene_list = [x.id for x in self.model.genes[:100]]
        
        results = single_deletion(self.model, gene_list)[0]
        fbc_results = single_deletion(fbc_model, gene_list)[0]
        _tolerance = 1e-6
        for gene in gene_list:
            _result = max(_tolerance, results[gene])
            _fbc_result = max(_tolerance, fbc_results[gene])
            self.assertAlmostEqual(_result, _fbc_result, places=3)

    
    @skipIf(not __test_matlab, "scipy.io.loadmat required")
    def test_mat_read_write(self):
        """read and write COBRA toolbox models in .mat files"""
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
