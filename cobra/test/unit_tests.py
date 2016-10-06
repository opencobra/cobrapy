import sys
from unittest import TestCase, TestLoader, TextTestRunner, skipIf
from copy import copy, deepcopy
from pickle import loads, dumps, HIGHEST_PROTOCOL
import warnings
import re

if __name__ == "__main__":
    sys.path.insert(0, "../..")
    from cobra.test import create_test_model
    from cobra import Object, Model, Metabolite, Reaction, DictList
    sys.path.pop(0)
else:
    from . import create_test_model
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
        # ensure trying to index with an object that is a different object
        # also raises an error
        self.assertRaises(ValueError, self.list.index, Object("test1"))

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
        self.assertIs(self.list.get_by_id("test2"), obj2)
        self.assertEqual(len(self.list), 2)

    def testInsert(self):
        obj2 = Object("a")
        self.list.insert(0, obj2)
        self.assertEqual(self.list.index(obj2), 0)
        self.assertEqual(self.list.index("test1"), 1)
        self.assertIs(self.list.get_by_id("a"), obj2)
        self.assertEqual(len(self.list), 2)
        self.assertRaises(ValueError, self.list.append, obj2)

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
        obj2.name = "foobar1"
        self.list.append(obj2)
        result = self.list.query("test1")  # matches only test1
        self.assertEqual(len(result), 1)
        self.assertEqual(result[0], self.obj)
        result = self.list.query("foo", "name")  # matches only test2
        self.assertEqual(len(result), 1)
        self.assertEqual(result[0], obj2)
        result = self.list.query("test")  # matches test1 and test2
        self.assertEqual(len(result), 2)
        # test with a regular expression
        result = self.list.query(re.compile("test[0-9]"))
        self.assertEqual(len(result), 2)
        result = self.list.query(re.compile("test[29]"))
        self.assertEqual(len(result), 1)
        # test query of name
        result = self.list.query(re.compile("foobar."), "name")
        self.assertEqual(len(result), 1)

    def testRemoval(self):
        obj_list = DictList(Object("test%d" % (i)) for i in range(2, 10))
        del obj_list[3]
        self.assertNotIn("test5", obj_list)
        self.assertEqual(obj_list.index(obj_list[-1]), len(obj_list) - 1)
        self.assertEqual(len(obj_list), 7)
        del obj_list[3:5]
        self.assertNotIn("test6", obj_list)
        self.assertNotIn("test7", obj_list)
        self.assertEqual(obj_list.index(obj_list[-1]), len(obj_list) - 1)
        self.assertEqual(len(obj_list), 5)
        removed = obj_list.pop(1)
        self.assertEqual(obj_list.index(obj_list[-1]), len(obj_list) - 1)
        self.assertEqual(removed.id, "test3")
        self.assertNotIn("test3", obj_list)
        self.assertEqual(len(obj_list), 4)
        removed = obj_list.pop()
        self.assertEqual(removed.id, "test9")
        self.assertNotIn(removed.id, obj_list)
        self.assertEqual(len(obj_list), 3)

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

    def testSortandReverse(self):
        dl = DictList(Object("test%d" % (i)) for i in reversed(range(10)))
        self.assertEqual(dl[0].id, "test9")
        dl.sort()
        self.assertEqual(len(dl), 10)
        self.assertEqual(dl[0].id, "test0")
        self.assertEqual(dl.index("test0"), 0)
        dl.reverse()
        self.assertEqual(dl[0].id, "test9")
        self.assertEqual(dl.index("test0"), 9)

    def testDir(self):
        """makes sure tab complete will work"""
        attrs = dir(self.list)
        self.assertIn("test1", attrs)
        self.assertIn("_dict", attrs)  # attribute of DictList

    def testUnion(self):
        self.list.union([Object("test1"), Object("test2")])
        # should only add 1 element
        self.assertEqual(len(self.list), 2)
        self.assertEqual(self.list.index("test2"), 1)


class CobraTestCase(TestCase):
    def setUp(self):
        self.model = create_test_model("textbook")
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
        # test ability to handle uppercase AND/OR
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            reaction.gene_reaction_rule = "(b1 AND b2) OR (b3 and b4)"
        self.assertEqual(reaction.gene_reaction_rule,
                         "(b1 and b2) or (b3 and b4)")
        self.assertEqual(len(reaction.genes), 4)
        # ensure regular expressions correctly extract genes from malformed
        # GPR string
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            reaction.gene_reaction_rule = "(a1 or a2"
            self.assertEqual(len(reaction.genes), 2)
            reaction.gene_reaction_rule = "(forT or "
            self.assertEqual(len(reaction.genes), 1)

    def testGPR_modification(self):
        model = self.model
        reaction = model.reactions.get_by_id("PGI")
        old_gene = list(reaction.genes)[0]
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
        fake_gene.name = "foo_gene"
        self.assertEqual(reaction.gene_name_reaction_rule, fake_gene.name)

    def test_add_metabolite(self):
        model = self.model
        reaction = model.reactions.get_by_id("PGI")
        reaction.add_metabolites({model.metabolites[0]: 1})
        self.assertIn(model.metabolites[0], reaction._metabolites)
        fake_metabolite = Metabolite("fake")
        reaction.add_metabolites({fake_metabolite: 1})
        self.assertIn(fake_metabolite, reaction._metabolites)
        self.assertTrue(model.metabolites.has_id("fake"))
        self.assertIs(model.metabolites.get_by_id("fake"), fake_metabolite)

        # test adding by string
        reaction.add_metabolites({"g6p_c": -1})  # already in reaction
        self.assertTrue(
            reaction._metabolites[model.metabolites.get_by_id("g6p_c")], -2)
        reaction.add_metabolites({"h_c": 1})  # not currently in reaction
        self.assertTrue(
            reaction._metabolites[model.metabolites.get_by_id("h_c")], 1)
        with self.assertRaises(KeyError):
            reaction.add_metabolites({"missing": 1})

        # test adding to a new Reaction
        reaction = Reaction("test")
        self.assertEqual(len(reaction._metabolites), 0)
        reaction.add_metabolites({Metabolite("test_met"): -1})
        self.assertEqual(len(reaction._metabolites), 1)

    def test_subtract_metabolite(self):
        model = self.model
        reaction = model.reactions.get_by_id("PGI")
        reaction.subtract_metabolites(reaction.metabolites)
        self.assertEqual(len(reaction.metabolites), 0)

    def test_mass_balance(self):
        model = self.model
        reaction = model.reactions.get_by_id("PGI")
        # should be balanced now
        self.assertEqual(len(reaction.check_mass_balance()), 0)
        # should not be balanced after adding a hydrogen
        reaction.add_metabolites({model.metabolites.get_by_id("h_c"): 1})
        imbalance = reaction.check_mass_balance()
        self.assertIn("H", imbalance)
        self.assertEqual(imbalance["H"], 1)

    def test_build_from_string(self):
        model = self.model
        m = len(model.metabolites)
        pgi = model.reactions.get_by_id("PGI")
        pgi.reaction = "g6p_c --> f6p_c"
        self.assertEqual(pgi.lower_bound, 0)
        pgi.bounds = (0, 1000)
        self.assertEqual(pgi.bounds, (0, 1000))
        self.assertEqual(pgi.reversibility, False)
        pgi.reaction = "g6p_c <== f6p_c"
        self.assertEqual(pgi.upper_bound, 0)
        self.assertEqual(pgi.reaction.strip(), "g6p_c <-- f6p_c")
        pgi.reaction = "g6p_c --> f6p_c + h2o_c"
        self.assertIn(model.metabolites.h2o_c, pgi._metabolites)
        pgi.build_reaction_from_string("g6p_c --> f6p_c + foo", verbose=False)
        self.assertNotIn(model.metabolites.h2o_c, pgi._metabolites)
        self.assertIn("foo", model.metabolites)
        self.assertIn(model.metabolites.foo, pgi._metabolites)
        self.assertEqual(len(model.metabolites), m + 1)

    def test_copy(self):
        model = self.model
        PGI = model.reactions.PGI
        copied = PGI.copy()
        self.assertIsNot(PGI, copied)
        self.assertIs(PGI._model, model)
        self.assertIsNot(copied._model, model)
        # the copy should refer to different metabolites and genes
        for met in copied.metabolites:
            self.assertIsNot(met, model.metabolites.get_by_id(met.id))
            self.assertIsNot(met.model, model)
        for gene in copied.genes:
            self.assertIsNot(gene, model.genes.get_by_id(gene.id))
            self.assertIsNot(gene.model, model)

    def test_iadd(self):
        model = self.model
        PGI = model.reactions.PGI
        EX_h2o = model.reactions.EX_h2o_e
        original_PGI_gpr = PGI.gene_reaction_rule
        PGI += EX_h2o
        self.assertEqual(PGI.gene_reaction_rule, original_PGI_gpr)
        self.assertEqual(PGI.metabolites[model.metabolites.h2o_e], -1.0)
        # original should not have changed
        self.assertEqual(EX_h2o.gene_reaction_rule, '')
        self.assertEqual(EX_h2o.metabolites[model.metabolites.h2o_e], -1.0)
        # what about adding a reaction not in the model
        new_reaction = Reaction("test")
        new_reaction.add_metabolites({Metabolite("A"): -1, Metabolite("B"): 1})
        PGI += new_reaction
        self.assertEqual(PGI.gene_reaction_rule, original_PGI_gpr)
        self.assertEqual(len(PGI.gene_reaction_rule), 5)
        # and vice versa
        new_reaction += PGI
        self.assertEqual(len(new_reaction.metabolites), 5)  # not 7
        self.assertEqual(len(new_reaction.genes), 1)
        self.assertEqual(new_reaction.gene_reaction_rule, original_PGI_gpr)
        # what about combining 2 gpr's
        model.reactions.ACKr += model.reactions.ACONTa
        self.assertEqual(model.reactions.ACKr.gene_reaction_rule,
                         "(b2296 or b3115 or b1849) and (b0118 or b1276)")
        self.assertEqual(len(model.reactions.ACKr.genes), 5)

    def test_add(self):
        # not in place addition should work on a copy
        model = self.model
        new = model.reactions.PGI + model.reactions.EX_h2o_e
        self.assertIsNot(new._model, model)
        self.assertEqual(len(new.metabolites), 3)
        # the copy should refer to different metabolites and genes
        # This currently fails because add_metabolites does not copy.
        # Should that be changed?
        # for met in new.metabolites:
        #    self.assertIsNot(met, model.metabolites.get_by_id(met.id))
        #    self.assertIsNot(met.model, model)
        for gene in new.genes:
            self.assertIsNot(gene, model.genes.get_by_id(gene.id))
            self.assertIsNot(gene.model, model)

    def test_mul(self):
        new = self.model.reactions.PGI * 2
        self.assertEqual(set(new.metabolites.values()), {-2, 2})

    def test_sub(self):
        model = self.model
        new = model.reactions.PGI - model.reactions.EX_h2o_e
        self.assertIsNot(new._model, model)
        self.assertEqual(len(new.metabolites), 3)


class TestCobraMetabolites(CobraTestCase):
    def test_metabolite_formula(self):
        met = Metabolite("water")
        met.formula = "H2O"
        self.assertEqual(met.elements, {"H": 2, "O": 1})
        self.assertEqual(met.formula_weight, 18.01528)

    def test_foruma_element_setting(self):
        model = self.model
        met = model.metabolites[1]
        orig_formula = str(met.formula)
        orig_elements = dict(met.elements)

        met.formula = ''
        self.assertEqual(met.elements, {})
        met.elements = orig_elements
        self.assertEqual(met.formula, orig_formula)


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
        # test adding a different metabolite with the same name as an
        # existing one uses the metabolite in the model
        r2 = Reaction("test_foo_reaction2")
        self.model.add_reaction(r2)
        r2.add_metabolites({Metabolite(self.model.metabolites[0].id): 1})
        self.assertIs(self.model.metabolites[0], list(r2._metabolites)[0])

    def test_add_reaction_from_other_model(self):
        model = self.model
        other = model.copy()
        for i in other.reactions:
            i.id += "_other"
        other.repair()
        model.add_reactions(other.reactions)
        # what if the other reaction has an error in its GPR
        m1 = create_test_model("textbook")
        m2 = create_test_model("textbook")
        m1.reactions.PGI.remove_from_model()
        m2.genes.b4025._reaction.clear()
        m1.add_reaction(m2.reactions.PGI)

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
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
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
        _model.add_reactions((x.copy() for x in self.model.reactions))
        _genes = []
        _metabolites = []
        for x in _model.reactions:
            _genes.extend(x.genes)
            _metabolites.extend(x._metabolites)

        orphan_genes = [x for x in _genes if x.model is not _model]
        orphan_metabolites = [x for x in _metabolites if x.model is not _model]
        self.assertEqual(len(orphan_genes), 0,
                         msg='It looks like there are dangling genes when '
                         'running Model.add_reactions')
        self.assertEqual(len(orphan_metabolites), 0,
                         msg='It looks like there are dangling metabolites '
                         'when running Model.add_reactions')

    def test_change_objective(self):
        biomass = self.model.reactions.get_by_id("Biomass_Ecoli_core")
        atpm = self.model.reactions.get_by_id("ATPM")
        self.model.objective = atpm.id
        self.assertEqual(atpm.objective_coefficient, 1.)
        self.assertEqual(biomass.objective_coefficient, 0.)
        self.assertEqual(self.model.objective, {atpm: 1.})
        # change it back using object itself
        self.model.objective = biomass
        self.assertEqual(atpm.objective_coefficient, 0.)
        self.assertEqual(biomass.objective_coefficient, 1.)
        # set both to 1 with a list
        self.model.objective = [atpm, biomass]
        self.assertEqual(atpm.objective_coefficient, 1.)
        self.assertEqual(biomass.objective_coefficient, 1.)
        # set both using a dict
        self.model.objective = {atpm: 0.2, biomass: 0.3}
        self.assertEqual(atpm.objective_coefficient, 0.2)
        self.assertEqual(biomass.objective_coefficient, 0.3)
        # test setting by index
        self.model.objective = self.model.reactions.index(atpm)
        self.assertEqual(self.model.objective, {atpm: 1.})
        # test by setting list of indexes
        self.model.objective = map(self.model.reactions.index, [atpm, biomass])
        self.assertEqual(self.model.objective, {atpm: 1., biomass: 1.})


@skipIf(scipy is None, "scipy required for ArrayBasedModel")
class TestCobraArrayModel(TestCobraModel):
    def setUp(self):
        model = create_test_model("textbook").to_array_based_model()
        self.model_class = model.__class__
        self.model = model

    def test_array_based_model(self):
        m = len(self.model.metabolites)
        n = len(self.model.reactions)
        assertEqual = self.assertEqual  # alias
        for matrix_type in ["scipy.dok_matrix", "scipy.lil_matrix"]:
            model = create_test_model("textbook").\
                to_array_based_model(matrix_type=matrix_type)
            assertEqual(model.S[7, 0], -1)
            assertEqual(model.S[43, 0], 0)
            model.S[43, 0] = 1
            assertEqual(model.S[43, 0], 1)
            assertEqual(
                model.reactions[0]._metabolites[model.metabolites[43]], 1)
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

            # test something for all the attributes
            model.lower_bounds[2] = -1
            assertEqual(model.reactions[2].lower_bound, -1)
            assertEqual(model.lower_bounds[2], -1)
            model.objective_coefficients[2] = 1
            assertEqual(model.reactions[2].objective_coefficient, 1)
            assertEqual(model.objective_coefficients[2], 1)
            model.b[2] = 1
            assertEqual(model.metabolites[2]._bound, 1)
            assertEqual(model.b[2], 1)
            model.constraint_sense[2] = "L"
            assertEqual(model.metabolites[2]._constraint_sense, "L")
            assertEqual(model.constraint_sense[2], "L")

            # test resize matrix on reaction removal
            m, n = model.S.shape
            model.remove_reactions([model.reactions[2]], remove_orphans=False)
            self.assertEqual(len(model.metabolites), model.S.shape[0])
            self.assertEqual(len(model.reactions), model.S.shape[1])
            self.assertEqual(model.S.shape, (m, n - 1))

    def test_array_based_model_add(self):
        m = len(self.model.metabolites)
        n = len(self.model.reactions)
        for matrix_type in ["scipy.dok_matrix", "scipy.lil_matrix"]:
            model = create_test_model("textbook").\
                to_array_based_model(matrix_type=matrix_type)
            test_reaction = Reaction("test")
            test_reaction.add_metabolites({model.metabolites[0]: 4})
            test_reaction.lower_bound = -3.14
            model.add_reaction(test_reaction)
            self.assertEqual(len(model.reactions), n + 1)
            self.assertEqual(model.S.shape, (m, n + 1))
            self.assertEqual(len(model.lower_bounds), n + 1)
            self.assertEqual(len(model.upper_bounds), n + 1)
            self.assertEqual(model.S[0, n], 4)
            self.assertEqual(model.S[7, 0], -1)
            self.assertEqual(model.lower_bounds[n], -3.14)

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

    def test_array_based_bounds_setting(self):
        model = self.model
        bounds = [0.0] * len(model.reactions)
        model.lower_bounds = bounds
        self.assertEqual(type(model.reactions[0].lower_bound), float)
        self.assertAlmostEqual(model.reactions[0].lower_bound, 0.0)
        model.upper_bounds[1] = 1234.0
        self.assertAlmostEqual(model.reactions[1].upper_bound, 1234.0)
        model.upper_bounds[9:11] = [100.0, 200.0]
        self.assertAlmostEqual(model.reactions[9].upper_bound, 100.0)
        self.assertAlmostEqual(model.reactions[10].upper_bound, 200.0)
        model.upper_bounds[9:11] = 123.0
        self.assertAlmostEqual(model.reactions[9].upper_bound, 123.0)
        self.assertAlmostEqual(model.reactions[10].upper_bound, 123.0)


# make a test suite to run all of the tests
loader = TestLoader()
suite = loader.loadTestsFromModule(sys.modules[__name__])


def test_all():
    TextTestRunner(verbosity=2).run(suite)

if __name__ == "__main__":
    test_all()
