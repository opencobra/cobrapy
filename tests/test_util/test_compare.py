"""Test functions of cobra.util.compare.py ."""

from cobra import Model
from cobra.util.compare import compare_reaction_state
import warnings
from typing import Iterable

import numpy as np
import pytest

from cobra import Gene
from cobra.core import GPR, Configuration, Metabolite, Model, Reaction

def test_reaction_copies_are_equivalent(model: Model) -> None:
    """Test that a copy of a reaction will return true when using compare functions."""
    reaction = model.reactions.get_by_id("PGI")
    old_reaction = reaction.copy()
    equivalent, comparison = compare_reaction_state(reaction, old_reaction)
    assert equivalent
    assert comparison['same'] == reaction.__getstate__().keys()
    assert comparison['modified'] == {}
    assert comparison['added'] == set()
    assert comparison['removed'] == set()

def test_reaction_different_ids(model: Model) -> None:
    """Test that reactions that differ in ids are not identical."""
    reaction = model.reactions.get_by_id("PGI")
    old_reaction = reaction.copy()
    reaction.id = 'PGI2'
    equivalent, comparison = compare_reaction_state(reaction, old_reaction)
    assert not equivalent
    assert comparison['same'] == set(reaction.__getstate__().keys()).difference({'_id'})
    assert comparison['modified'] == {'_id': ('PGI2', 'PGI')}
    assert comparison['added'] == set()
    assert comparison['removed'] == set()


def test_reaction_different_names(model: Model) -> None:
    """Test that reactions that differ in names are not identical."""
    reaction = model.reactions.get_by_id("PGI")
    old_reaction = reaction.copy()
    reaction.name = reaction.name + ' '
    equivalent, comparison = compare_reaction_state(reaction, old_reaction)
    assert not equivalent
    assert comparison['same'] == set(reaction.__getstate__().keys()).difference({'name'})
    assert comparison['modified'] == {'name': (old_reaction.name + ' ', old_reaction.name)}
    assert comparison['added'] == set()
    assert comparison['removed'] == set()
    reaction.name = None
    equivalent, comparison = compare_reaction_state(reaction, old_reaction)
    assert not equivalent
    assert comparison['same'] == set(reaction.__getstate__().keys()).difference({'name'})
    assert comparison['modified'] == {'name': (None, old_reaction.name)}
    assert comparison['added'] == set()
    assert comparison['removed'] == set()
    reaction.name = ''
    equivalent, comparison = compare_reaction_state(reaction, old_reaction)
    assert not equivalent
    assert comparison['same'] == set(reaction.__getstate__().keys()).difference({'name'})
    assert comparison['modified'] == {'name': ('', old_reaction.name)}
    assert comparison['added'] == set()
    assert comparison['removed'] == set()
    reaction.name = 'Test'
    equivalent, comparison = compare_reaction_state(reaction, old_reaction)
    assert not equivalent
    assert comparison['same'] == set(reaction.__getstate__().keys()).difference({'name'})
    assert comparison['modified'] == {'name': ('Test', old_reaction.name)}
    assert comparison['added'] == set()
    assert comparison['removed'] == set()


def test_reaction_gpr_modification(model: Model) -> None:
    """Test reactions are not equivalent after GPR/gene rule manipulations."""
    reaction = model.reactions.get_by_id("PGI")
    old_reaction = reaction.copy()
    old_gene = list(reaction.genes)[0]

    # Add an existing 'gene' to the GPR
    reaction.gene_reaction_rule = "s0001"
    equivalent, comparison = compare_reaction_state(reaction, old_reaction)
    assert not equivalent
    assert comparison['modified']['_genes'] == ({'s0001'}, {old_gene.id})
    assert comparison['modified']['_gpr'] == (reaction.gpr, old_reaction.gpr)
    assert comparison['same'] == set(reaction.__getstate__().keys()).difference({'_gpr', '_genes'})

    old_reaction = reaction.copy()
    old_reaction.gpr = GPR()
    equivalent, comparison = compare_reaction_state(reaction, old_reaction)
    assert not equivalent
    assert comparison['modified']['_genes'] == ({'s0001'}, set())
    assert comparison['modified']['_gpr'] == (reaction.gpr, GPR())
    assert comparison['same'] == set(reaction.__getstate__().keys()).difference(
        {'_gpr', '_genes'})



def test_compare_rxn_bounds(model: Model) -> None:
    """Test reaction bounds setting for a scenario."""
    acald_reaction = model.reactions.ACALD
    new_reaction = acald_reaction.copy()
    new_reaction.bounds = (
        acald_reaction.lower_bound - 100,
        acald_reaction.lower_bound - 100,
    )
    equivalent, comparison = compare_reaction_state(acald_reaction, new_reaction)
    assert not equivalent
    assert comparison['same'] == set(acald_reaction.__getstate__().keys()).difference(
        {'_lower_bound', '_upper_bound'})
    assert comparison['modified'] == {'_lower_bound': (-1000, -1100),
                                      '_upper_bound': (1000, -1100)}
    assert comparison['added'] == set()
    assert comparison['removed'] == set()
    new_reaction = acald_reaction.copy()
    new_reaction.bounds = (
        acald_reaction.upper_bound + 100,
        acald_reaction.upper_bound + 100,
    )
    equivalent, comparison = compare_reaction_state(acald_reaction, new_reaction)
    assert not equivalent
    assert comparison['same'] == set(acald_reaction.__getstate__().keys()).difference(
        {'_lower_bound', '_upper_bound'})
    assert comparison['modified'] == {'_lower_bound': (-1000, 1100),
                                      '_upper_bound': (1000, 1100)}
    assert comparison['added'] == set()
    assert comparison['removed'] == set()
    new_reaction = acald_reaction.copy()
    new_reaction.lower_bound = -100
    equivalent, comparison = compare_reaction_state(acald_reaction, new_reaction)
    assert not equivalent
    assert comparison['same'] == set(acald_reaction.__getstate__().keys()).difference(
        {'_lower_bound'})
    assert comparison['modified'] == {'_lower_bound': (-1000, -100)}
    assert comparison['added'] == set()
    assert comparison['removed'] == set()
    new_reaction = acald_reaction.copy()
    new_reaction.upper_bound = 100
    equivalent, comparison = compare_reaction_state(acald_reaction, new_reaction)
    assert not equivalent
    assert comparison['same'] == set(acald_reaction.__getstate__().keys()).difference(
        {'_upper_bound'})
    assert comparison['modified'] == {'_upper_bound': (1000, 100)}
    assert comparison['added'] == set()
    assert comparison['removed'] == set()
    new_reaction = acald_reaction.copy()
    new_reaction.knock_out()
    equivalent, comparison = compare_reaction_state(acald_reaction, new_reaction)
    assert not equivalent
    assert comparison['same'] == set(acald_reaction.__getstate__().keys()).difference(
        {'_lower_bound', '_upper_bound'})
    assert comparison['modified'] == {'_lower_bound': (-1000, 0),
                                      '_upper_bound': (1000, 0)}
    assert comparison['added'] == set()
    assert comparison['removed'] == set()
    new_reaction = acald_reaction.copy()
    new_reaction.bounds = (0, 0)
    equivalent, comparison = compare_reaction_state(acald_reaction, new_reaction)
    assert not equivalent
    assert comparison['same'] == set(acald_reaction.__getstate__().keys()).difference(
        {'_lower_bound', '_upper_bound'})
    assert comparison['modified'] == {'_lower_bound': (-1000, 0),
                                      '_upper_bound': (1000, 0)}
    assert comparison['added'] == set()
    assert comparison['removed'] == set()


def test_compare_reaction_different_subsystem(model: Model) -> None:
    """Test that differences in subsystem are picked up by compare function."""

# {'_id', '_genes', '_gpr', '_lower_bound', '_annotation', 'subsystem', '_metabolites', '_model', 'name', 'notes', '_upper_bound'}


def test_add_metabolite_from_solved_model(solved_model: Model) -> None:
    """Test metabolite addition to a reaction from a solved model."""
    solution, model = solved_model
    pgi_reaction = model.reactions.PGI
    test_met = model.metabolites[0]
    pgi_reaction.add_metabolites({test_met: 42}, combine=False)
    assert pgi_reaction.metabolites[test_met] == 42.0
    assert (
        model.constraints[test_met.id].expression.as_coefficients_dict()[
            pgi_reaction.forward_variable
        ]
        == 42.0
    )
    assert (
        model.constraints[test_met.id].expression.as_coefficients_dict()[
            pgi_reaction.reverse_variable
        ]
        == -42.0
    )

    pgi_reaction.add_metabolites({test_met: -10}, combine=True)
    assert pgi_reaction.metabolites[test_met] == 32.0
    assert (
        model.constraints[test_met.id].expression.as_coefficients_dict()[
            pgi_reaction.forward_variable
        ]
        == 32.0
    )
    assert (
        model.constraints[test_met.id].expression.as_coefficients_dict()[
            pgi_reaction.reverse_variable
        ]
        == -32.0
    )

    pgi_reaction.add_metabolites({test_met: 0}, combine=False)
    with pytest.raises(KeyError):
        pgi_reaction.metabolites[test_met]
    assert (
        model.constraints[test_met.id].expression.as_coefficients_dict()[
            pgi_reaction.forward_variable
        ]
        == 0
    )
    assert (
        model.constraints[test_met.id].expression.as_coefficients_dict()[
            pgi_reaction.reverse_variable
        ]
        == 0
    )


def test_add_metabolite(model: Model) -> None:
    """Test metabolite addition to a reaction from an unsolved model."""
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
        assert reaction._metabolites[model.metabolites.get_by_id("g6p_c")] == -2
        reaction.add_metabolites({"h_c": 1})
        assert reaction._metabolites[model.metabolites.get_by_id("h_c")] == 1
        with pytest.raises(KeyError):
            reaction.add_metabolites({"missing": 1})

    assert reaction._metabolites[model.metabolites.get_by_id("g6p_c")] == -1
    assert model.metabolites.h_c not in reaction._metabolites

    # Test combine=False
    reaction = model.reactions.get_by_id("ATPM")
    old_stoich = reaction._metabolites[model.metabolites.get_by_id("h2o_c")]
    with model:
        reaction.add_metabolites({"h2o_c": 2.5}, combine=False)
        assert reaction._metabolites[model.metabolites.get_by_id("h2o_c")] == 2.5

    assert reaction._metabolites[model.metabolites.get_by_id("h2o_c")] == old_stoich

    # Test adding to a new Reaction
    reaction = Reaction("test")
    assert len(reaction._metabolites) == 0
    reaction.add_metabolites({Metabolite("test_met"): -1})
    assert len(reaction._metabolites) == 1


def test_subtract_metabolite(model: Model, solver: Iterable) -> None:
    """Test metabolite deletion from a reaction associated with an unsolved model."""
    reaction = model.reactions.get_by_id("PGI")
    reaction.subtract_metabolites(reaction.metabolites)


def test_mass_balance(model: Model) -> None:
    """Test mass balance of metabolites of a reaction."""
    reaction = model.reactions.get_by_id("PGI")
    # Should be balanced now
    assert len(reaction.check_mass_balance()) == 0
    # Should not be balanced after adding a hydrogen
    reaction.add_metabolites({model.metabolites.get_by_id("h_c"): 1})
    imbalance = reaction.check_mass_balance()
    assert "H" in imbalance
    assert imbalance["H"] == 1


def test_iadd(model: Model) -> None:
    """Test in-place addition of reaction."""
    PGI = model.reactions.PGI
    EX_h2o = model.reactions.EX_h2o_e
    original_PGI_gpr = PGI.gene_reaction_rule
    PGI += EX_h2o
    assert PGI.gene_reaction_rule == original_PGI_gpr
    assert PGI.metabolites[model.metabolites.h2o_e] == -1.0
    # Original should not change
    assert EX_h2o.gene_reaction_rule == ""
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
    expected_rule = "(b2296 or b3115 or b1849) and (b0118 or b1276)"
    assert model.reactions.ACKr.gene_reaction_rule == expected_rule
    assert len(model.reactions.ACKr.genes) == 5


def test_iadd_with_context(model: Model) -> None:
    """Test in-place addition of reaction is reversed with context."""
    PGI = model.reactions.PGI
    EX_h2o = model.reactions.EX_h2o_e
    original_PGI_gene_reaction_rule = PGI.gene_reaction_rule
    with model:
        PGI += EX_h2o
        assert PGI.gene_reaction_rule == original_PGI_gene_reaction_rule
        assert PGI.metabolites[model.metabolites.h2o_e] == -1.0
    assert PGI.gene_reaction_rule == original_PGI_gene_reaction_rule
    assert model.metabolites.h2o_e not in PGI.metabolites.keys()
    # Add a reaction not in the model
    new_reaction = Reaction("test")
    new_reaction.add_metabolites({Metabolite("A"): -1, Metabolite("B"): 1})
    with model:
        PGI += new_reaction
    assert PGI.gene_reaction_rule == original_PGI_gene_reaction_rule
    assert len(PGI.gene_reaction_rule) == 5
    # Combine two GPRs
    expected_rule = "(b2296 or b3115 or b1849) and (b0118 or b1276)"
    old_rule = model.reactions.ACKr.gene_reaction_rule
    with model:
        model.reactions.ACKr += model.reactions.ACONTa
        assert model.reactions.ACKr.gene_reaction_rule == expected_rule
        assert len(model.reactions.ACKr.genes) == 5
    assert model.reactions.ACKr.gene_reaction_rule == old_rule
    assert old_rule != expected_rule
    assert len(model.reactions.ACKr.genes) == 3


def test_radd(model: Model) -> None:
    """Test __radd__ for a reaction."""
    new = sum([model.reactions.PGI, model.reactions.EX_h2o_e])
    assert new._model is not model
    assert len(new.metabolites) == 3


def test_mul(model: Model) -> None:
    """Test scalar multiplication of factors with a reaction."""
    new = model.reactions.PGI * 2
    assert set(new.metabolites.values()) == {-2, 2}


def test_sub(model: Model) -> None:
    """Test reaction subtraction."""
    new = model.reactions.PGI - model.reactions.EX_h2o_e
    assert new._model is not model
    assert len(new.metabolites) == 3

def test_add_metabolites_combine_true(model: Model) -> None:
    """Test metabolite addition to reaction (with combine = True)."""
    test_metabolite = Metabolite("test")
    for reaction in model.reactions:
        reaction.add_metabolites({test_metabolite: -66}, combine=True)
        assert reaction.metabolites[test_metabolite] == -66
        assert (
            model.constraints["test"].get_linear_coefficients(
                [reaction.forward_variable]
            )[reaction.forward_variable]
            == -66
        )
        assert (
            model.constraints["test"].get_linear_coefficients(
                [reaction.reverse_variable]
            )[reaction.reverse_variable]
            == 66
        )
        already_included_metabolite = list(reaction.metabolites.keys())[0]
        previous_coefficient = reaction.get_coefficient(already_included_metabolite.id)
        reaction.add_metabolites({already_included_metabolite: 10}, combine=True)
        new_coefficient = previous_coefficient + 10
        assert reaction.metabolites[already_included_metabolite] == new_coefficient
        assert (
            model.constraints[already_included_metabolite.id].get_linear_coefficients(
                [reaction.forward_variable]
            )[reaction.forward_variable]
            == new_coefficient
        )
        assert (
            model.constraints[already_included_metabolite.id].get_linear_coefficients(
                [reaction.reverse_variable]
            )[reaction.reverse_variable]
            == -new_coefficient
        )

def test_reaction_imul(model: Model) -> None:
    """Test in-place scalar factor multiplication to reaction."""
    with model:
        model.reactions.EX_glc__D_e *= 100
        assert (
            model.constraints.glc__D_e.expression.coeff(model.variables.EX_glc__D_e)
            == -100.0
        )
        assert model.reactions.EX_glc__D_e.reaction == "100.0 glc__D_e <=> "

    assert (
        model.constraints.glc__D_e.expression.coeff(model.variables.EX_glc__D_e) == -1.0
    )
    assert model.reactions.EX_glc__D_e.reaction == "glc__D_e <=> "

    with model:
        model.reactions.EX_glc__D_e *= -2
        assert model.reactions.EX_glc__D_e.bounds == (-1000.0, 10.0)
        assert model.reactions.EX_glc__D_e.reaction == " <=> 2.0 glc__D_e"

    assert model.reactions.EX_glc__D_e.bounds == (-10, 1000.0)
    assert model.reactions.EX_glc__D_e.reaction == "glc__D_e <=> "


def test_compartment_changes(model: Model) -> None:
    """Test reaction compartment change."""
    rxn = model.reactions.EX_ac_e
    assert rxn.reactants[0].compartment in rxn.compartments
    rxn.reactants[0].compartment = "blub"
    assert rxn.reactants[0].compartment in rxn.compartments

## Test model
def test_add(model: Model) -> None:
    """Test reaction addition to model."""
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


def test_removal_from_model_retains_bounds(model: Model) -> None:
    """Test reaction removal from a model, retains its bounds."""
    model_cp = model.copy()
    reaction = model_cp.reactions.ACALD
    assert reaction.model == model_cp
    assert reaction.lower_bound == -1000.0
    assert reaction.upper_bound == 1000.0
    assert reaction._lower_bound == -1000.0
    assert reaction._upper_bound == 1000.0
    model_cp.remove_reactions([reaction])
    assert reaction.model is None
    assert reaction.lower_bound == -1000.0
    assert reaction.upper_bound == 1000.0
    assert reaction._lower_bound == -1000.0
    assert reaction._upper_bound == 1000.0


def test_remove_from_model(model: Model) -> None:
    """Test reaction removal from model."""
    pgi = model.reactions.PGI
    g6p = model.metabolites.g6p_c
    pgi_flux = model.optimize().fluxes["PGI"]
    assert abs(pgi_flux) > 1e-6

    with model:
        pgi.remove_from_model()
        assert pgi.model is None
        assert "PGI" not in model.reactions
        assert pgi.id not in model.variables
        assert pgi.reverse_id not in model.variables
        assert pgi not in g6p.reactions
        model.optimize()

    assert "PGI" in model.reactions
    assert pgi.id in model.variables
    assert pgi.reverse_id in model.variables
    assert pgi.forward_variable.problem is model.solver
    assert pgi in g6p.reactions
    assert g6p in pgi.metabolites
    assert np.isclose(pgi_flux, model.optimize().fluxes["PGI"])

