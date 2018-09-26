# -*- coding: utf-8 -*-

"""Test functions of reaction.py"""

from __future__ import absolute_import

import warnings

import pytest

from cobra.core import Metabolite, Model, Reaction

stable_optlang = ["glpk", "cplex", "gurobi"]


def test_gpr():
    model = Model()
    reaction = Reaction("test")

    # Set GPR to a reaction not in a model
    reaction.gene_reaction_rule = "(g1 or g2) and g3"
    assert reaction.gene_reaction_rule == "(g1 or g2) and g3"
    assert len(reaction.genes) == 3

    # Adding reaction with a GPR propagates to the model
    model.add_reactions([reaction])
    assert len(model.genes) == 3

    # Ensure the gene objects are the same in the model and reaction
    reaction_gene = list(reaction.genes)[0]
    model_gene = model.genes.get_by_id(reaction_gene.id)
    assert reaction_gene is model_gene

    # Test ability to handle uppercase AND/OR
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        reaction.gene_reaction_rule = "(b1 AND b2) OR (b3 and b4)"
        assert reaction.gene_reaction_rule == "(b1 and b2) or (b3 and b4)"
        assert len(reaction.genes) == 4

    # Ensure regular expressions correctly extract genes from malformed
    # GPR string
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        reaction.gene_reaction_rule = "(a1 or a2"
        assert len(reaction.genes) == 2
        reaction.gene_reaction_rule = "(forT or "
        assert len(reaction.genes) == 1


def test_gpr_modification(model):
    reaction = model.reactions.get_by_id("PGI")
    old_gene = list(reaction.genes)[0]
    new_gene = model.genes.get_by_id("s0001")

    # Add an existing 'gene' to the GPR
    reaction.gene_reaction_rule = 's0001'
    assert new_gene in reaction.genes
    assert reaction in new_gene.reactions

    # Remove old gene correctly
    assert old_gene not in reaction.genes
    assert reaction not in old_gene.reactions

    # Add a new 'gene' to the GPR
    reaction.gene_reaction_rule = 'fake_gene'
    assert model.genes.has_id("fake_gene")
    fake_gene = model.genes.get_by_id("fake_gene")
    assert fake_gene in reaction.genes
    assert reaction in fake_gene.reactions
    fake_gene.name = "foo_gene"
    assert reaction.gene_name_reaction_rule == fake_gene.name


def test_gene_knock_out(model):
    rxn = Reaction('rxn')
    rxn.add_metabolites({Metabolite('A'): -1, Metabolite('B'): 1})
    rxn.gene_reaction_rule = 'A2B1 or A2B2 and A2B3'
    assert hasattr(list(rxn.genes)[0], 'knock_out')
    model.add_reaction(rxn)
    with model:
        model.genes.A2B1.knock_out()
        assert not model.genes.A2B1.functional
        model.genes.A2B3.knock_out()
        assert not rxn.functional
    assert model.genes.A2B3.functional
    assert rxn.functional
    model.genes.A2B1.knock_out()
    assert not model.genes.A2B1.functional
    assert model.reactions.rxn.functional
    model.genes.A2B3.knock_out()
    assert not model.reactions.rxn.functional


def test_str():
    rxn = Reaction('rxn')
    rxn.add_metabolites({Metabolite('A'): -1, Metabolite('B'): 1})
    assert str(rxn) == 'rxn: A --> B'


@pytest.mark.parametrize("solver", stable_optlang)
def test_add_metabolite_benchmark(model, benchmark, solver):
    reaction = model.reactions.get_by_id("PGI")
    many_metabolites = dict((m, 1) for m in model.metabolites[0:50])

    def add_remove_metabolite():
        reaction.add_metabolites(many_metabolites)
        if not getattr(model, 'solver', None):
            solver_dict[solver].create_problem(model)
        for m, c in many_metabolites.items():
            try:
                reaction.subtract_metabolites(
                    {m: reaction.get_coefficient(m)})
            except KeyError:
                pass

    benchmark(add_remove_metabolite)


def test_add_metabolite(model):
    with pytest.raises(ValueError):
        model.add_metabolites(Metabolite())
    with model:
        with model:
            reaction = model.reactions.get_by_id("PGI")
            reaction.add_metabolites({model.metabolites[0]: 1})
            assert model.metabolites[0] in reaction._metabolites
            fake_metabolite = Metabolite("fake")
            reaction.add_metabolites({fake_metabolite: 1})
            assert fake_metabolite in reaction._metabolites
            assert model.metabolites.has_id("fake")
            assert model.metabolites.get_by_id("fake") is fake_metabolite
            assert len(model._contexts[0]._history) == 0

    assert fake_metabolite._model is None
    assert fake_metabolite not in reaction._metabolites
    assert "fake" not in model.metabolites

    # Test adding by string
    with model:
        reaction.add_metabolites({"g6p_c": -1})  # already in reaction
        assert reaction._metabolites[
            model.metabolites.get_by_id("g6p_c")] == -2
        reaction.add_metabolites({"h_c": 1})
        assert reaction._metabolites[
            model.metabolites.get_by_id("h_c")] == 1
        with pytest.raises(KeyError):
            reaction.add_metabolites({"missing": 1})

    assert reaction._metabolites[
        model.metabolites.get_by_id("g6p_c")] == -1
    assert model.metabolites.h_c not in reaction._metabolites

    # Test combine=False
    reaction = model.reactions.get_by_id("ATPM")
    old_stoich = reaction._metabolites[
        model.metabolites.get_by_id("h2o_c")]
    with model:
        reaction.add_metabolites({'h2o_c': 2.5}, combine=False)
        assert reaction._metabolites[
            model.metabolites.get_by_id("h2o_c")] == 2.5

    assert reaction._metabolites[
        model.metabolites.get_by_id("h2o_c")] == old_stoich

    # Test adding to a new Reaction
    reaction = Reaction("test")
    assert len(reaction._metabolites) == 0
    reaction.add_metabolites({Metabolite("test_met"): -1})
    assert len(reaction._metabolites) == 1


@pytest.mark.parametrize("solver", stable_optlang)
def test_subtract_metabolite_benchmark(model, benchmark, solver):
    benchmark(test_subtract_metabolite, model, solver)


@pytest.mark.parametrize("solver", stable_optlang)
def test_subtract_metabolite(model, solver):
    reaction = model.reactions.get_by_id("PGI")
    reaction.subtract_metabolites(reaction.metabolites)
    if not getattr(model, 'solver', None):
        solver_dict[solver].create_problem(model)
        assert len(reaction.metabolites) == 0


def test_mass_balance(model):
    reaction = model.reactions.get_by_id("PGI")
    # Should be balanced now
    assert len(reaction.check_mass_balance()) == 0
    # Should not be balanced after adding a hydrogen
    reaction.add_metabolites({model.metabolites.get_by_id("h_c"): 1})
    imbalance = reaction.check_mass_balance()
    assert "H" in imbalance
    assert imbalance["H"] == 1


def test_build_from_string(model):
    m = len(model.metabolites)
    pgi = model.reactions.get_by_id("PGI")
    old_bounds = pgi.bounds

    with model:
        pgi.reaction = "g6p_c --> f6p_c"
        assert pgi.lower_bound == 0

    assert pgi.bounds == old_bounds

    pgi.bounds = (0, 1000)
    assert pgi.bounds == (0, 1000)
    assert not pgi.reversibility
    pgi.reaction = "g6p_c <== f6p_c"
    assert pgi.upper_bound == 0
    assert pgi.reaction.strip() == "g6p_c <-- f6p_c"
    pgi.reaction = "g6p_c --> f6p_c + h2o_c"
    assert model.metabolites.h2o_c, pgi._metabolites

    with model:
        pgi.build_reaction_from_string("g6p_c --> f6p_c + foo",
                                       verbose=False)
        assert model.metabolites.h2o_c not in pgi._metabolites
        assert "foo" in model.metabolites
        assert model.metabolites.foo in pgi._metabolites
        assert len(model.metabolites) == m + 1

    assert model.metabolites.h2o_c in pgi._metabolites
    assert "foo" not in model.metabolites
    with pytest.raises(AttributeError):
        assert model.metabolites.foo
        assert len(model.metabolites) == m


def test_bounds_setter(model):
    rxn = model.reactions.get_by_id("PGI")
    with pytest.raises(AssertionError):
        rxn.bounds = (1, 0)


def test_copy(model):
    PGI = model.reactions.PGI
    copied = PGI.copy()
    assert PGI is not copied
    assert PGI._model is model
    assert copied._model is not model
    # The copy should refer to different metabolites and genes
    for met in copied.metabolites:
        assert met is not model.metabolites.get_by_id(met.id)
        assert met.model is not model
    for gene in copied.genes:
        assert gene is not model.genes.get_by_id(gene.id)
        assert gene.model is not model


def test_iadd(model):
    PGI = model.reactions.PGI
    EX_h2o = model.reactions.EX_h2o_e
    original_PGI_gpr = PGI.gene_reaction_rule
    PGI += EX_h2o
    assert PGI.gene_reaction_rule == original_PGI_gpr
    assert PGI.metabolites[model.metabolites.h2o_e] == -1.0
    # Original should not change
    assert EX_h2o.gene_reaction_rule == ''
    assert EX_h2o.metabolites[model.metabolites.h2o_e] == -1.0
    # Add a reaction not in the model
    new_reaction = Reaction("test")
    new_reaction.add_metabolites({Metabolite("A"): -1, Metabolite("B"): 1})
    PGI += new_reaction
    assert PGI.gene_reaction_rule == original_PGI_gpr
    assert len(PGI.gene_reaction_rule) == 5
    # And vice versa
    new_reaction += PGI
    assert len(new_reaction.metabolites) == 5  # not
    assert len(new_reaction.genes) == 1
    assert new_reaction.gene_reaction_rule == original_PGI_gpr
    # Combine two GPRs
    model.reactions.ACKr += model.reactions.ACONTa
    expected_rule = '(b2296 or b3115 or b1849) and (b0118 or b1276)'
    assert model.reactions.ACKr.gene_reaction_rule == expected_rule
    assert len(model.reactions.ACKr.genes) == 5


def test_add(model):
    # Not in place addition should work on a copy
    new = model.reactions.PGI + model.reactions.EX_h2o_e
    assert new._model is not model
    assert len(new.metabolites) == 3
    # The copy should refer to different metabolites and genes
    # This currently fails because add_metabolites does not copy.
    # Should that be changed?
    # for met in new.metabolites:
    #    assert met is not model.metabolites.get_by_id(met.id)
    #    assert met.model is not model
    for gene in new.genes:
        assert gene is not model.genes.get_by_id(gene.id)
        assert gene.model is not model


def test_radd(model):
    new = sum([model.reactions.PGI, model.reactions.EX_h2o_e])
    assert new._model is not model
    assert len(new.metabolites) == 3


def test_mul(model):
    new = model.reactions.PGI * 2
    assert set(new.metabolites.values()) == {-2, 2}


def test_sub(model):
    new = model.reactions.PGI - model.reactions.EX_h2o_e
    assert new._model is not model
    assert len(new.metabolites) == 3


def test_repr_html_(model):
    assert '<table>' in model.reactions[0]._repr_html_()
