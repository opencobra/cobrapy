import pytest

from cobra.core import Metabolite, Model, Reaction
from cobra.manipulation import *
from .conftest import model, salmonella


class TestManipulation:
    """Test functions in cobra.manipulation"""

    def test_canonical_form(self, model):
        solver = 'cglpk'
        # add G constraint to test
        g_constr = Metabolite("SUCCt2_2__test_G_constraint")
        g_constr._constraint_sense = "G"
        g_constr._bound = 5.0
        model.reactions.get_by_id("SUCCt2_2").add_metabolites({g_constr: 1})
        assert abs(model.optimize("maximize", solver=solver).f - 0.855) < 0.001
        # convert to canonical form
        model = canonical_form(model)
        assert abs(
            model.optimize("maximize", solver=solver).f - 0.855) < 10 ** -3

    def test_canonical_form_minimize(self, model):
        solver = 'cglpk'
        # make a minimization problem
        model.reactions.get_by_id("Biomass_Ecoli_core").lower_bound = 0.5
        for reaction in model.reactions:
            reaction.objective_coefficient = reaction.id == "GAPD"
        assert abs(
            model.optimize("minimize", solver=solver).f - 6.27) < 10 ** -3
        # convert to canonical form. Convert minimize to maximize
        model = canonical_form(model, objective_sense="minimize")
        assert abs(
            model.optimize("maximize", solver=solver).f + 6.27) < 10 ** -3
        # lower bounds should now be <= constraints
        assert model.reactions.get_by_id(
            "Biomass_Ecoli_core").lower_bound == 0.0

    def test_modify_reversible(self, model):
        solver = 'cglpk'
        model1 = model.copy()
        model1.optimize(solver=solver)
        model2 = model.copy()
        convert_to_irreversible(model2)
        model2.optimize(solver=solver)
        assert abs(model1.solution.f - model2.solution.f) < 10 ** -3
        revert_to_reversible(model2)
        model2.optimize(solver=solver)
        assert abs(model1.solution.f - model2.solution.f) < 10 ** -3
        # Ensure revert_to_reversible is robust to solutions generated both
        # before and after reversibility conversion, or not solved at all.
        model3 = model.copy()
        model3.optimize(solver=solver)
        convert_to_irreversible(model3)
        revert_to_reversible(model3)
        assert abs(model1.solution.f - model3.solution.f) < 10 ** -3
        # test reaction where both bounds are negative
        model4 = model.copy()
        glc = model4.reactions.get_by_id("EX_glc__D_e")
        glc.upper_bound = -1
        convert_to_irreversible(model4)
        model4.optimize(solver=solver)
        assert abs(model1.solution.f - model4.solution.f) < 10 ** -3
        glc_rev = model4.reactions.get_by_id(glc.notes["reflection"])
        assert glc_rev.lower_bound == 1
        assert glc.upper_bound == 0
        revert_to_reversible(model4)
        assert glc.upper_bound == -1

    def test_escape_ids(self, model):
        model.reactions.PGI.gene_reaction_rule = "a.b or c"
        assert "a.b" in model.genes
        escape_ID(model)
        assert "a.b" not in model.genes

    def test_rename_gene(self, model):
        original_name = model.genes.b1241.name
        rename_dict = {"b1241": "foo", "hello": "world",
                       "b2465": "b3919", "bar": "2935"}
        modify.rename_genes(model, rename_dict)
        for i in rename_dict:
            assert i not in model.genes
        assert "foo" in model.genes
        # make sure the object name was preserved
        assert model.genes.foo.name == original_name
        # make sure the reactions are correct
        assert len(model.genes.foo.reactions) == 2
        assert model.reactions.ACALD.gene_reaction_rule == "b0351 or foo"
        assert model.reactions.TPI.gene_reaction_rule == "b3919"
        assert model.reactions.TPI.genes == {model.genes.b3919}
        assert model.reactions.TKT1.gene_reaction_rule == "b2935 or b3919"
        assert model.reactions.TKT1.genes == {
            model.genes.b2935, model.genes.b3919}
        assert model.genes.b3919.reactions == {model.reactions.get_by_id(i)
                                               for i in
                                               ("TKT1", "TKT2", "TPI")}

    def test_gene_knockout_computation(self, salmonella):
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
            assert removed1 == expected_reactions
            assert removed2 == expected_reactions
            delete_model_genes(m, gene_ids, cumulative_deletions=False)
            assert get_removed(m) == expected_reaction_ids
            undelete_model_genes(m)

        gene_list = ['STM1067', 'STM0227']
        dependent_reactions = {'3HAD121', '3HAD160', '3HAD80', '3HAD140',
                               '3HAD180', '3HAD100', '3HAD181', '3HAD120',
                               '3HAD60', '3HAD141', '3HAD161', 'T2DECAI',
                               '3HAD40'}
        test_computation(salmonella, gene_list, dependent_reactions)
        test_computation(salmonella, ['STM4221'], {'PGI'})
        test_computation(salmonella, ['STM1746.S'], {'4PEPTabcpp'})
        # test cumulative behavior
        delete_model_genes(salmonella, gene_list[:1])
        delete_model_genes(salmonella, gene_list[1:],
                           cumulative_deletions=True)
        delete_model_genes(salmonella, ["STM4221"],
                           cumulative_deletions=True)
        dependent_reactions.add('PGI')
        assert get_removed(salmonella) == dependent_reactions
        # non-cumulative following cumulative
        delete_model_genes(salmonella, ["STM4221"],
                           cumulative_deletions=False)
        assert get_removed(salmonella) == {'PGI'}
        # make sure on reset that the bounds are correct
        reset_bound = salmonella.reactions.get_by_id("T2DECAI").upper_bound
        assert reset_bound == 1000.
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
        assert len(m.reactions) == 8
        rxns = m.reactions
        rxns.r1.gene_reaction_rule = "(a and b) or (c and a)"
        rxns.r2.gene_reaction_rule = "(a and b and d and e)"
        rxns.r3.gene_reaction_rule = "(a and b) or (b and c)"
        rxns.r4.gene_reaction_rule = "(f and b) or (b and c)"
        rxns.r5.gene_reaction_rule = "x"
        rxns.r6.gene_reaction_rule = "y"
        rxns.r7.gene_reaction_rule = "x or     z"
        rxns.r8.gene_reaction_rule = ""
        assert "a" in m.genes
        assert "x" in m.genes
        remove_genes(m, ["a"], remove_reactions=False)
        assert "a" not in m.genes
        assert "x" in m.genes
        assert rxns.r1.gene_reaction_rule == ""
        assert rxns.r2.gene_reaction_rule == ""
        assert rxns.r3.gene_reaction_rule == "b and c"
        assert rxns.r4.gene_reaction_rule == "(f and b) or (b and c)"
        assert rxns.r5.gene_reaction_rule == "x"
        assert rxns.r6.gene_reaction_rule == "y"
        assert rxns.r7.genes == {m.genes.x, m.genes.z}
        assert rxns.r8.gene_reaction_rule == ""
        remove_genes(m, ["x"], remove_reactions=True)
        assert len(m.reactions) == 7
        assert "r5" not in m.reactions
        assert "x" not in m.genes
        assert rxns.r1.gene_reaction_rule == ""
        assert rxns.r2.gene_reaction_rule == ""
        assert rxns.r3.gene_reaction_rule == "b and c"
        assert rxns.r4.gene_reaction_rule == "(f and b) or (b and c)"
        assert rxns.r6.gene_reaction_rule == "y"
        assert rxns.r7.gene_reaction_rule == "z"
        assert rxns.r7.genes == {m.genes.z}
        assert rxns.r8.gene_reaction_rule == ""

    def test_sbo_annotation(self, model):
        rxns = model.reactions
        rxns.EX_o2_e.annotation.clear()
        fake_DM = Reaction("DM_h_c")
        model.add_reaction(fake_DM)
        fake_DM.add_metabolites({model.metabolites.get_by_id("h_c"): -1})
        # this exchange will be set wrong. The function should not overwrite
        # an existing SBO annotation
        rxns.get_by_id("EX_h_e").annotation["SBO"] = "SBO:0000628"
        add_SBO(model)
        assert rxns.EX_o2_e.annotation["SBO"] == "SBO:0000627"
        assert rxns.DM_h_c.annotation["SBO"] == "SBO:0000628"
        assert rxns.EX_h_e.annotation["SBO"] == "SBO:0000628"

    def test_validate_formula_compartment(self, model):
        model.metabolites[1].compartment = "fake"
        model.metabolites[1].formula = "(a*.bcde)"
        errors = check_metabolite_compartment_formula(model)
        assert len(errors) == 2

    def test_validate_mass_balance(self, model):
        assert len(check_mass_balance(model)) == 0
        # if we remove the SBO term which marks the reaction as
        # mass balanced, then the reaction should be detected as
        # no longer mass balanced
        EX_rxn = model.reactions.query(lambda r: r.boundary)[0]
        EX_rxn.annotation.pop("SBO")
        balance = check_mass_balance(model)
        assert len(balance) == 1
        assert EX_rxn in balance
        m1 = Metabolite('m1', formula='()')
        r1 = Reaction('r1')
        r1.add_metabolites({m1: 1})
        with pytest.raises(ValueError):
            r1.check_mass_balance()

    def test_prune_unused(self, model):
        metabolite = model.metabolites.ru5p__D_c
        [model.reactions.get_by_id(x).remove_from_model() for x in
         ['RPI', 'RPE', 'GND']]
        unused = delete.prune_unused_metabolites(model)
        assert unused[0] is metabolite
        reaction = Reaction('foo')
        model.add_reaction(reaction)
        unused = delete.prune_unused_reactions(model)
        assert unused[0] is reaction
