import sys
from unittest import TestCase, TestLoader, TextTestRunner, skipIf
from copy import copy, deepcopy
from pickle import loads, dumps, HIGHEST_PROTOCOL

if __name__ == "__main__":
    sys.path.insert(0, "../..")
    from cobra.test import data_directory, create_test_model
    from cobra.test import ecoli_mat, ecoli_pickle
    from cobra.test import salmonella_sbml, salmonella_pickle
    from cobra import Object, Model, Metabolite, Reaction, DictList
    sys.path.pop(0)
else:
    from . import data_directory, create_test_model
    from . import ecoli_mat, ecoli_pickle
    from . import salmonella_sbml, salmonella_pickle
    from .. import Object, Model, Metabolite, Reaction, DictList

# libraries which may or may not be installed
libraries = ["scipy"]
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

    def testContains(self):
        self.assertIn(self.obj, self.list)
        self.assertIn(self.obj.id, self.list)
        self.assertNotIn(Object("not_in"), self.list)
        self.assertNotIn('not_in', self.list)

    def testIndex(self):
        self.assertEqual(self.list.index("test1"), 0)
        self.assertEqual(self.list.index(self.obj), 0)
        self.assertRaises(ValueError, self.list.index, "f")
        self.assertRaises(ValueError, self.list.index, Object("f"))

    def testIndependent(self):
        a = DictList([Object("o1"), Object("o2")])
        b = DictList()
        self.assertIn("o1", a)
        self.assertNotIn("o1", b)
        b.append(Object("o3"))
        self.assertNotIn("o3", a)
        self.assertIn("o3", b)

    def testAppend(self):
        obj2 = Object("test2")
        self.list.append(obj2)
        self.assertRaises(ValueError, self.list.append, Object("test1"))
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
        self.assertRaises(ValueError, self.list.extend, [Object("test1")])
        # Even if the object is unique, if it is present twice in the new
        # list, it should still raise an exception
        self.assertRaises(ValueError, self.list.extend,
                          [Object("testd"), Object("testd")])

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
        self.assertIsNot(sum, self.list)
        self.assertIsNot(sum, obj_list)
        self.assertEqual(self.list[0].id, "test1")
        self.assertEqual(sum[1].id, "test2")
        self.assertEqual(sum.get_by_id("test2"), obj_list[0])
        self.assertEqual(sum[8].id, "test9")
        self.assertEqual(len(self.list), 1)
        self.assertEqual(len(sum), 9)

    def testInitCopy(self):
        self.list.append(Object("test2"))
        copied = DictList(self.list)
        self.assertIsNot(self.list, copied)
        self.assertIsInstance(copied, self.list.__class__)
        self.assertEqual(len(self.list), len(copied))
        for i, v in enumerate(self.list):
            self.assertEqual(self.list[i].id, copied[i].id)
            self.assertEqual(i, copied.index(v.id))
            self.assertIs(self.list[i], copied[i])
            self.assertIs(v, copied.get_by_id(v.id))

    def testSlice(self):
        self.list.append(Object("test2"))
        self.list.append(Object("test3"))
        sliced = self.list[:-1]
        self.assertIsNot(self.list, sliced)
        self.assertIsInstance(sliced, self.list.__class__)
        self.assertEqual(len(self.list), len(sliced) + 1)
        for i, v in enumerate(sliced):
            self.assertEqual(self.list[i].id, sliced[i].id)
            self.assertEqual(i, sliced.index(v.id))
            self.assertIs(self.list[i], sliced[i])
            self.assertIs(self.list[i], sliced.get_by_id(v.id))

    def testCopy(self):
        self.list.append(Object("test2"))
        copied = copy(self.list)
        self.assertIsNot(self.list, copied)
        self.assertIsInstance(copied, self.list.__class__)
        self.assertEqual(len(self.list), len(copied))
        for i, v in enumerate(self.list):
            self.assertEqual(self.list[i].id, copied[i].id)
            self.assertEqual(i, copied.index(v.id))
            self.assertIs(self.list[i], copied[i])
            self.assertIs(v, copied.get_by_id(v.id))

    def testDeepcopy(self):
        self.list.append(Object("test2"))
        copied = deepcopy(self.list)
        self.assertIsNot(self.list, copied)
        self.assertIsInstance(copied, self.list.__class__)
        self.assertEqual(len(self.list), len(copied))
        for i, v in enumerate(self.list):
            self.assertEqual(self.list[i].id, copied[i].id)
            self.assertEqual(i, copied.index(v.id))
            self.assertIsNot(self.list[i], copied[i])
            self.assertIsNot(v, copied.get_by_id(v.id))

    def testPickle(self):
        self.list.append(Object("test2"))
        for protocol in range(HIGHEST_PROTOCOL):
            pickle_str = dumps(self.list, protocol=protocol)
            copied = loads(pickle_str)
            self.assertIsNot(self.list, copied)
            self.assertIsInstance(copied, self.list.__class__)
            self.assertEqual(len(self.list), len(copied))
            for i, v in enumerate(self.list):
                self.assertEqual(self.list[i].id, copied[i].id)
                self.assertEqual(i, copied.index(v.id))
                self.assertIsNot(self.list[i], copied[i])
                self.assertIsNot(v, copied.get_by_id(v.id))

    def testQuery(self):
        obj2 = Object("test2")
        self.list.append(obj2)
        result = self.list.query("test1")  # matches only test1
        self.assertEqual(len(result), 1)
        self.assertEqual(result[0], self.obj)
        result = self.list.query("test")  # matches test1 and test2
        self.assertEqual(len(result), 2)

    def testRemoval(self):
        obj_list = DictList(Object("test%d" % (i)) for i in range(2, 10))
        del obj_list[3]
        self.assertNotIn("test5", obj_list)
        self.assertEqual(obj_list.index(obj_list[-1]), len(obj_list) - 1)
        del obj_list[3:5]
        self.assertNotIn("test6", obj_list)
        self.assertNotIn("test7", obj_list)
        self.assertEqual(obj_list.index(obj_list[-1]), len(obj_list) - 1)
        removed = obj_list.pop(1)
        self.assertEqual(obj_list.index(obj_list[-1]), len(obj_list) - 1)
        self.assertEqual(removed.id, "test3")
        self.assertNotIn("test3", obj_list)

    def testSet(self):
        obj_list = DictList(Object("test%d" % (i)) for i in range(10))
        obj_list[4] = Object("testa")
        self.assertEqual(obj_list.index("testa"), 4)
        self.assertEqual(obj_list[4].id, "testa")
        obj_list[5:7] = [Object("testb"), Object("testc")]
        self.assertEqual(obj_list.index("testb"), 5)
        self.assertEqual(obj_list[5].id, "testb")
        self.assertEqual(obj_list.index("testc"), 6)
        self.assertEqual(obj_list[6].id, "testc")
        # Even if the object is unique, if it is present twice in the new
        # list, it should still raise an exception
        self.assertRaises(ValueError, obj_list.__setitem__, slice(5, 7),
                          [Object("testd"), Object("testd")])


class CobraTestCase(TestCase):
    def setUp(self):
        self.model = create_test_model()
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

    def test_build_from_string(self):
        model = self.model
        pgi = model.reactions.get_by_id("PGI")
        pgi.reaction = "g6p_c --> f6p_c"
        self.assertEqual(pgi.lower_bound, 0)
        pgi.reaction = "g6p_c <== f6p_c"
        self.assertEqual(pgi.upper_bound, 0)
        pgi.reaction = "g6p_c --> f6p_c + h2o_c"
        self.assertIn(model.metabolites.h2o_c, pgi._metabolites)
        pgi.build_reaction_from_string("g6p_c --> f6p_c + foo", verbose=False)
        self.assertNotIn(model.metabolites.h2o_c, pgi._metabolites)
        self.assertIn("foo", model.metabolites)
        self.assertIn(model.metabolites.foo, pgi._metabolites)
        self.assertEqual(len(model.metabolites), 1803)


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
        copy_in_model = self.model.metabolites.get_by_id(copy_metabolite.id)
        self.assertIsNot(copy_metabolite, copy_in_model)
        self.assertIs(type(copy_in_model), Metabolite)
        self.assertTrue(dummy_reaction in actual_metabolite._reaction)

    def test_add_reaction_from_other_model(self):
        model = self.model
        other = model.copy()
        for i in other.reactions:
            i.id += "_other"
        other.reactions._generate_index()
        model.add_reactions(other.reactions)

    def test_model_remove_reaction(self):
        old_reaction_count = len(self.model.reactions)
        self.model.remove_reactions(["PGI"])
        self.assertEqual(len(self.model.reactions), old_reaction_count - 1)
        with self.assertRaises(KeyError):
            self.model.reactions.get_by_id("PGI")
        self.model.remove_reactions(self.model.reactions[:1])
        self.assertEqual(len(self.model.reactions), old_reaction_count - 2)
        tmp_metabolite = Metabolite("testing")
        self.model.reactions[0].add_metabolites({tmp_metabolite: 1})
        self.assertIn(tmp_metabolite, self.model.metabolites)
        self.model.remove_reactions(self.model.reactions[:1],
                                    remove_orphans=True)
        self.assertNotIn(tmp_metabolite, self.model.metabolites)

    def test_reaction_remove(self):
        model = self.model
        old_reaction_count = len(model.reactions)
        tmp_metabolite = Metabolite("testing")
        # Delete without removing orphan
        model.reactions[0].add_metabolites({tmp_metabolite: 1})
        self.assertEqual(len(tmp_metabolite.reactions), 1)
        # esnsure the stoichiometry is still the same using different objects
        removed_reaction = model.reactions[0]
        original_stoich = {i.id: value for i, value
                           in removed_reaction._metabolites.items()}
        model.reactions[0].remove_from_model(remove_orphans=False)
        self.assertEqual(len(original_stoich),
                         len(removed_reaction._metabolites))
        for met in removed_reaction._metabolites:
            self.assertEqual(original_stoich[met.id],
                             removed_reaction._metabolites[met])
            self.assertIsNot(met, model.metabolites)
        # make sure it's still in the model
        self.assertIn(tmp_metabolite, model.metabolites)
        self.assertEqual(len(tmp_metabolite.reactions), 0)
        self.assertEqual(len(self.model.reactions), old_reaction_count - 1)

        # Now try it with removing orphans
        model.reactions[0].add_metabolites({tmp_metabolite: 1})
        self.assertEqual(len(tmp_metabolite.reactions), 1)
        model.reactions[0].remove_from_model(remove_orphans=True)
        self.assertNotIn(tmp_metabolite, model.metabolites)
        self.assertEqual(len(tmp_metabolite.reactions), 0)
        self.assertEqual(len(self.model.reactions), old_reaction_count - 2)

        # It shouldn't remove orphans if it's in 2 reactions however
        model.reactions[0].add_metabolites({tmp_metabolite: 1})
        model.reactions[1].add_metabolites({tmp_metabolite: 1})
        self.assertEqual(len(tmp_metabolite.reactions), 2)
        model.reactions[0].remove_from_model(remove_orphans=False)
        self.assertIn(tmp_metabolite, model.metabolites)
        self.assertEqual(len(tmp_metabolite.reactions), 1)
        self.assertEqual(len(self.model.reactions), old_reaction_count - 3)

    def test_reaction_delete(self):
        model = self.model
        old_reaction_count = len(model.reactions)
        tmp_metabolite = Metabolite("testing")
        # Delete without removing orphan
        model.reactions[0].add_metabolites({tmp_metabolite: 1})
        self.assertEqual(len(tmp_metabolite.reactions), 1)
        model.reactions[0].delete(remove_orphans=False)
        # make sure it's still in the model
        self.assertIn(tmp_metabolite, model.metabolites)
        self.assertEqual(len(tmp_metabolite.reactions), 0)
        self.assertEqual(len(self.model.reactions), old_reaction_count - 1)

        # Now try it with removing orphans
        model.reactions[0].add_metabolites({tmp_metabolite: 1})
        self.assertEqual(len(tmp_metabolite.reactions), 1)
        model.reactions[0].delete(remove_orphans=True)
        self.assertNotIn(tmp_metabolite, model.metabolites)
        self.assertEqual(len(tmp_metabolite.reactions), 0)
        self.assertEqual(len(self.model.reactions), old_reaction_count - 2)

        # It shouldn't remove orphans if it's in 2 reactions however
        model.reactions[0].add_metabolites({tmp_metabolite: 1})
        model.reactions[1].add_metabolites({tmp_metabolite: 1})
        self.assertEqual(len(tmp_metabolite.reactions), 2)
        model.reactions[0].delete(remove_orphans=False)
        self.assertIn(tmp_metabolite, model.metabolites)
        self.assertEqual(len(tmp_metabolite.reactions), 1)
        self.assertEqual(len(self.model.reactions), old_reaction_count - 3)

    def test_remove_gene(self):
        target_gene = self.model.genes[0]
        gene_reactions = list(target_gene.reactions)
        target_gene.remove_from_model()
        self.assertEqual(target_gene.model, None)
        # make sure the reaction was removed from the model
        self.assertNotIn(target_gene, self.model.genes)
        # ensure the old reactions no longer have a record of the gene
        for reaction in gene_reactions:
            self.assertNotIn(target_gene, reaction.genes)

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
        for reaction, reaction_copy in zip(self.model.reactions,
                                           model_copy.reactions):
            self.assertEqual(reaction.id, reaction_copy.id)
            metabolites = sorted(i.id for i in reaction._metabolites)
            metabolites_copy = sorted(i.id for i in reaction_copy._metabolites)
            self.assertEqual(metabolites, metabolites_copy)

    def test_add_reaction_orphans(self):
        """test reaction addition

        Need to verify that no orphan genes or metabolites are
        contained in reactions after adding them to the model.
        """
        _model = self.model_class('test')
        _model.add_reactions([x.copy() for x in self.model.reactions])
        _genes = []
        _metabolites = []
        for x in _model.reactions:
            _genes.extend(x.genes)
            _metabolites.extend(x.metabolites)

        orphan_genes = [x for x in _genes if x.model is not _model]
        orphan_metabolites = [x for x in _metabolites if x.model is not _model]
        self.assertEqual(len(orphan_genes), 0,
                         msg='It looks like there are dangling genes when '
                         'running Model.add_reactions')
        self.assertEqual(len(orphan_metabolites), 0,
                         msg='It looks like there are dangling metabolites '
                         'when running Model.add_reactions')

    def test_change_objective(self):
        biomass = self.model.reactions.get_by_id("biomass_iRR1083_metals")
        atpm = self.model.reactions.get_by_id("ATPM")
        self.model.change_objective(atpm.id)
        self.assertEqual(atpm.objective_coefficient, 1.)
        self.assertEqual(biomass.objective_coefficient, 0.)
        # change it back using object itself
        self.model.change_objective(biomass)
        self.assertEqual(atpm.objective_coefficient, 0.)
        self.assertEqual(biomass.objective_coefficient, 1.)
        # set both to 1 with a list
        self.model.change_objective([atpm, biomass])
        self.assertEqual(atpm.objective_coefficient, 1.)
        self.assertEqual(biomass.objective_coefficient, 1.)
        # set both using a dict
        self.model.change_objective({atpm: 0.2, biomass: 0.3})
        self.assertEqual(atpm.objective_coefficient, 0.2)
        self.assertEqual(biomass.objective_coefficient, 0.3)


@skipIf(scipy is None, "scipy required for ArrayBasedModel")
class TestCobraArrayModel(TestCobraModel):
    def setUp(self):
        model = create_test_model().to_array_based_model()
        self.model_class = model.__class__
        self.model = model

    def test_array_based_model(self):
        assertEqual = self.assertEqual  # alias
        for matrix_type in ["scipy.dok_matrix", "scipy.lil_matrix"]:
            model = create_test_model().\
                to_array_based_model(matrix_type=matrix_type)
            self.assertEqual(model.S[1605, 0], -1)
            self.assertEqual(model.S[43, 0], 0)
            model.S[43, 0] = 1
            self.assertEqual(model.S[43, 0], 1)
            self.assertEqual(
                model.reactions[0].metabolites[model.metabolites[43]], 1)
            model.S[43, 0] = 0
            assertEqual(model.lower_bounds[0], model.reactions[0].lower_bound)
            assertEqual(model.lower_bounds[5], model.reactions[5].lower_bound)
            assertEqual(model.upper_bounds[0], model.reactions[0].upper_bound)
            assertEqual(model.upper_bounds[5], model.reactions[5].upper_bound)
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
            model = create_test_model().\
                to_array_based_model(matrix_type=matrix_type)
            test_reaction = Reaction("test")
            test_reaction.add_metabolites({model.metabolites[0]: 4})
            test_reaction.lower_bound = -3.14
            model.add_reaction(test_reaction)
            self.assertEqual(len(model.reactions), 2547)
            self.assertEqual(model.S.shape[1], 2547)
            self.assertEqual(len(model.lower_bounds), 2547)
            self.assertEqual(model.S[0, 2546], 4)
            self.assertEqual(model.S[1605, 0], -1)
            self.assertEqual(model.lower_bounds[2546], -3.14)

    def test_array_based_select(self):
        model = self.model
        atpm_select = model.reactions[model.lower_bounds > 0]
        self.assertEqual(len(atpm_select), 1)
        self.assertEqual(atpm_select[0].id, "ATPM")
        self.assertEqual(len(model.reactions[model.lower_bounds <= 0]),
                         len(model.reactions) - 1)
        # mismatched dimensions should give an error
        with self.assertRaises(TypeError):
            model.reactions[[True, False]]


# make a test suite to run all of the tests
loader = TestLoader()
suite = loader.loadTestsFromModule(sys.modules[__name__])


def test_all():
    TextTestRunner(verbosity=2).run(suite)

if __name__ == "__main__":
    test_all()
