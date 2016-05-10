from unittest import TestCase, TestLoader, TextTestRunner

import sys

if __name__ == "__main__":
    sys.path.insert(0, "../..")
    from cobra.test import create_test_model, data_directory
    from cobra.core import Metabolite, Model, Reaction
    from cobra.manipulation import *
    sys.path.pop(0)
else:
    from . import create_test_model, data_directory
    from ..core import Metabolite, Model, Reaction
    from ..manipulation import *


class TestManipulation(TestCase):
    """Test functions in cobra.manipulation"""

    def test_canonical_form(self):
        model = create_test_model("textbook")
        # add G constraint to test
        g_constr = Metabolite("SUCCt2_2__test_G_constraint")
        g_constr._constraint_sense = "G"
        g_constr._bound = 5.0
        model.reactions.get_by_id("SUCCt2_2").add_metabolites({g_constr: 1})
        self.assertAlmostEqual(model.optimize("maximize").f, 0.855, places=3)
        # convert to canonical form
        model = canonical_form(model)
        self.assertAlmostEqual(model.optimize("maximize").f, 0.855, places=3)

    def test_canonical_form_minimize(self):
        model = create_test_model("textbook")
        # make a minimization problem
        model.reactions.get_by_id("Biomass_Ecoli_core").lower_bound = 0.5
        for reaction in model.reactions:
            reaction.objective_coefficient = reaction.id == "GAPD"
        self.assertAlmostEqual(model.optimize("minimize").f, 6.27, places=3)
        # convert to canonical form. Convert minimize to maximize
        model = canonical_form(model, objective_sense="minimize")
        self.assertAlmostEqual(model.optimize("maximize").f, -6.27, places=3)
        # lower bounds should now be <= constraints
        self.assertEqual(
            model.reactions.get_by_id("Biomass_Ecoli_core").lower_bound, 0.0)

    def test_modify_reversible(self):
        model1 = create_test_model("textbook")
        model1.optimize()
        model2 = create_test_model("textbook")
        convert_to_irreversible(model2)
        model2.optimize()
        self.assertAlmostEqual(model1.solution.f, model2.solution.f, places=3)
        revert_to_reversible(model2)
        model2.optimize()
        self.assertAlmostEqual(model1.solution.f, model2.solution.f, places=3)

        # Ensure revert_to_reversible is robust to solutions generated both
        # before and after reversibility conversion, or not solved at all.
        model3 = create_test_model("textbook")
        model3.optimize()
        convert_to_irreversible(model3)
        revert_to_reversible(model3)
        self.assertAlmostEqual(model1.solution.f, model3.solution.f, places=3)

        # test reaction where both bounds are negative
        model4 = create_test_model("textbook")
        glc = model4.reactions.get_by_id("EX_glc__D_e")
        glc.upper_bound = -1
        convert_to_irreversible(model4)
        model4.optimize()
        self.assertAlmostEqual(model1.solution.f, model4.solution.f, places=3)
        glc_rev = model4.reactions.get_by_id(glc.notes["reflection"])
        self.assertEqual(glc_rev.lower_bound, 1)
        self.assertEqual(glc.upper_bound, 0)
        revert_to_reversible(model4)
        self.assertEqual(glc.upper_bound, -1)

    def test_escape_ids(self):
        model = create_test_model('textbook')
        model.reactions.PGI.gene_reaction_rule = "a.b or c"
        self.assertIn("a.b", model.genes)
        escape_ID(model)
        self.assertNotIn("a.b", model.genes)

    def test_rename_gene(self):
        model = create_test_model('textbook')
        original_name = model.genes.b1241.name
        rename_dict = {"b1241": "foo", "hello": "world",
                       "b2465": "b3919", "bar": "2935"}
        modify.rename_genes(model, rename_dict)
        for i in rename_dict:
            self.assertNotIn(i, model.genes)
        self.assertIn("foo", model.genes)
        # make sure the object name was preserved
        self.assertEqual(model.genes.foo.name, original_name)
        # make sure the reactions are correct
        self.assertEqual(len(model.genes.foo.reactions), 2)
        self.assertEqual(model.reactions.ACALD.gene_reaction_rule,
                         "b0351 or foo")
        self.assertEqual(model.reactions.TPI.gene_reaction_rule, "b3919")
        self.assertEqual(model.reactions.TPI.genes, {model.genes.b3919})
        self.assertEqual(model.reactions.TKT1.gene_reaction_rule,
                         "b2935 or b3919")
        self.assertEqual(model.reactions.TKT1.genes,
                         {model.genes.b2935, model.genes.b3919})
        self.assertEqual(model.genes.b3919.reactions,
                         {model.reactions.get_by_id(i)
                          for i in ("TKT1", "TKT2", "TPI")})

    def test_gene_knockout_computation(self):
        cobra_model = create_test_model()

        def find_gene_knockout_reactions_fast(cobra_model, gene_list):
            compiled_rules = get_compiled_gene_reaction_rules(
                cobra_model)
            return find_gene_knockout_reactions(
                cobra_model, gene_list,
                compiled_gene_reaction_rules=compiled_rules)

        def get_removed(m):
            return {x.id for x in m._trimmed_reactions}

        def test_computation(m, gene_ids, expected_reaction_ids):
            genes = [m.genes.get_by_id(i) for i in gene_ids]
            expected_reactions = {m.reactions.get_by_id(i)
                                  for i in expected_reaction_ids}
            removed1 = set(find_gene_knockout_reactions(m, genes))
            removed2 = set(find_gene_knockout_reactions_fast(m, genes))
            self.assertEqual(removed1, expected_reactions)
            self.assertEqual(removed2, expected_reactions)
            delete_model_genes(m, gene_ids, cumulative_deletions=False)
            self.assertEqual(get_removed(m), expected_reaction_ids)
            undelete_model_genes(m)

        gene_list = ['STM1067', 'STM0227']
        dependent_reactions = {'3HAD121', '3HAD160', '3HAD80', '3HAD140',
                               '3HAD180', '3HAD100', '3HAD181', '3HAD120',
                               '3HAD60', '3HAD141', '3HAD161', 'T2DECAI',
                               '3HAD40'}
        test_computation(cobra_model, gene_list, dependent_reactions)
        test_computation(cobra_model, ['STM4221'], {'PGI'})
        test_computation(cobra_model, ['STM1746.S'], {'4PEPTabcpp'})
        # test cumulative behavior
        delete_model_genes(cobra_model, gene_list[:1])
        delete_model_genes(cobra_model, gene_list[1:],
                           cumulative_deletions=True)
        delete_model_genes(cobra_model, ["STM4221"],
                           cumulative_deletions=True)
        dependent_reactions.add('PGI')
        self.assertEqual(get_removed(cobra_model), dependent_reactions)
        # non-cumulative following cumulative
        delete_model_genes(cobra_model, ["STM4221"],
                           cumulative_deletions=False)
        self.assertEqual(get_removed(cobra_model), {'PGI'})
        # make sure on reset that the bounds are correct
        reset_bound = cobra_model.reactions.get_by_id("T2DECAI").upper_bound
        self.assertEqual(reset_bound, 1000.)
        # test computation when gene name is a subset of another
        test_model = Model()
        test_reaction_1 = Reaction("test1")
        test_reaction_1.gene_reaction_rule = "eggs or (spam and eggspam)"
        test_model.add_reaction(test_reaction_1)
        test_computation(test_model, ["eggs"], set())
        test_computation(test_model, ["eggs", "spam"], {'test1'})
        # test computation with nested boolean expression
        test_reaction_1.gene_reaction_rule = \
            "g1 and g2 and (g3 or g4 or (g5 and g6))"
        test_computation(test_model, ["g3"], set())
        test_computation(test_model, ["g1"], {'test1'})
        test_computation(test_model, ["g5"], set())
        test_computation(test_model, ["g3", "g4", "g5"], {'test1'})
        # test computation when gene names are python expressions
        test_reaction_1.gene_reaction_rule = "g1 and (for or in)"
        test_computation(test_model, ["for", "in"], {'test1'})
        test_computation(test_model, ["for"], set())
        test_reaction_1.gene_reaction_rule = "g1 and g2 and g2.conjugate"
        test_computation(test_model, ["g2"], {"test1"})
        test_computation(test_model, ["g2.conjugate"], {"test1"})
        test_reaction_1.gene_reaction_rule = "g1 and (try:' or 'except:1)"
        test_computation(test_model, ["try:'"], set())
        test_computation(test_model, ["try:'", "'except:1"], {"test1"})

    def test_remove_genes(self):
        m = Model("test")
        m.add_reactions([Reaction("r" + str(i + 1)) for i in range(8)])
        self.assertEqual(len(m.reactions), 8)
        rxns = m.reactions
        rxns.r1.gene_reaction_rule = "(a and b) or (c and a)"
        rxns.r2.gene_reaction_rule = "(a and b and d and e)"
        rxns.r3.gene_reaction_rule = "(a and b) or (b and c)"
        rxns.r4.gene_reaction_rule = "(f and b) or (b and c)"
        rxns.r5.gene_reaction_rule = "x"
        rxns.r6.gene_reaction_rule = "y"
        rxns.r7.gene_reaction_rule = "x or     z"
        rxns.r8.gene_reaction_rule = ""
        self.assertIn("a", m.genes)
        self.assertIn("x", m.genes)
        remove_genes(m, ["a"], remove_reactions=False)
        self.assertNotIn("a", m.genes)
        self.assertIn("x", m.genes)
        self.assertEqual(rxns.r1.gene_reaction_rule, "")
        self.assertEqual(rxns.r2.gene_reaction_rule, "")
        self.assertEqual(rxns.r3.gene_reaction_rule, "b and c")
        self.assertEqual(rxns.r4.gene_reaction_rule, "(f and b) or (b and c)")
        self.assertEqual(rxns.r5.gene_reaction_rule, "x")
        self.assertEqual(rxns.r6.gene_reaction_rule, "y")
        self.assertEqual(rxns.r7.genes, {m.genes.x, m.genes.z})
        self.assertEqual(rxns.r8.gene_reaction_rule, "")
        remove_genes(m, ["x"], remove_reactions=True)
        self.assertEqual(len(m.reactions), 7)
        self.assertNotIn("r5", m.reactions)
        self.assertNotIn("x", m.genes)
        self.assertEqual(rxns.r1.gene_reaction_rule, "")
        self.assertEqual(rxns.r2.gene_reaction_rule, "")
        self.assertEqual(rxns.r3.gene_reaction_rule, "b and c")
        self.assertEqual(rxns.r4.gene_reaction_rule, "(f and b) or (b and c)")
        self.assertEqual(rxns.r6.gene_reaction_rule, "y")
        self.assertEqual(rxns.r7.gene_reaction_rule, "z")
        self.assertEqual(rxns.r7.genes, {m.genes.z})
        self.assertEqual(rxns.r8.gene_reaction_rule, "")

    def test_SBO_annotation(self):
        model = create_test_model("textbook")
        rxns = model.reactions
        rxns.EX_o2_e.annotation.clear()
        fake_DM = Reaction("DM_h_c")
        model.add_reaction(fake_DM)
        fake_DM.add_metabolites({model.metabolites.get_by_id("h_c"): -1})
        # this exchange will be set wrong. The function should not overwrite
        # an existing SBO annotation
        rxns.get_by_id("EX_h_e").annotation["SBO"] = "SBO:0000628"
        add_SBO(model)
        self.assertEqual(rxns.EX_o2_e.annotation["SBO"], "SBO:0000627")
        self.assertEqual(rxns.DM_h_c.annotation["SBO"], "SBO:0000628")
        self.assertEqual(rxns.EX_h_e.annotation["SBO"], "SBO:0000628")

    def test_validate_reaction_bounds(self):
        model = create_test_model("textbook")
        model.reactions[0].lower_bound = float("-inf")
        model.reactions[1].lower_bound = float("nan")
        model.reactions[0].upper_bound = float("inf")
        model.reactions[1].upper_bound = float("nan")
        errors = check_reaction_bounds(model)
        self.assertEqual(len(errors), 4)

    def test_validate_formula_compartment(self):
        model = create_test_model("textbook")
        model.metabolites[1].compartment = "fake"
        model.metabolites[1].formula = "(a*.bcde)"
        errors = check_metabolite_compartment_formula(model)
        self.assertEqual(len(errors), 2)

    def test_validate_mass_balance(self):
        model = create_test_model("textbook")
        self.assertEqual(len(check_mass_balance(model)), 0)
        # if we remove the SBO term which marks the reaction as
        # mass balanced, then the reaction should be detected as
        # no longer mass balanced
        EX_rxn = model.reactions.query("EX")[0]
        EX_rxn.annotation.pop("SBO")
        balance = check_mass_balance(model)
        self.assertEqual(len(balance), 1)
        self.assertIn(EX_rxn, balance)


# make a test suite to run all of the tests
loader = TestLoader()
suite = loader.loadTestsFromModule(sys.modules[__name__])


def test_all():
    TextTestRunner(verbosity=2).run(suite)

if __name__ == "__main__":
    test_all()
