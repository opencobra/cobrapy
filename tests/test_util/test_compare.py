"""Test Comparing functions in cobra.util.compare.py ."""


import numpy as np
import pytest

from cobra.core import GPR, Metabolite, Model, Reaction
from cobra.util.compare import compare_reaction_state, compare_state


def test_reaction_copies_are_equivalent(model: Model) -> None:
    """Test that a copy of a reaction will return true when using compare functions."""
    reaction = model.reactions.get_by_id("PGI")
    old_reaction = reaction.copy()
    equivalent, comparison = compare_reaction_state(reaction, old_reaction)
    assert equivalent
    assert comparison["same"] == reaction.__getstate__().keys()
    assert comparison["modified"] == {}
    assert comparison["added"] == set()
    assert comparison["removed"] == set()


def test_reaction_with_added_field_is_different(model: Model) -> None:
    """Test that reaction with an added field will be identified by compare."""
    reaction = model.reactions.get_by_id("PGI")
    old_reaction = reaction.copy()
    reaction.blah = None
    equivalent, comparison = compare_reaction_state(reaction, old_reaction)
    assert not equivalent
    assert comparison["same"] == old_reaction.__getstate__().keys()
    assert comparison["modified"] == {}
    assert comparison["added"] == {"blah"}
    assert comparison["removed"] == set()
    equivalent, comparison = compare_reaction_state(old_reaction, reaction)
    assert not equivalent
    assert comparison["same"] == old_reaction.__getstate__().keys()
    assert comparison["modified"] == {}
    assert comparison["added"] == set()
    assert comparison["removed"] == {"blah"}


def test_reaction_different_ids(model: Model) -> None:
    """Test that reactions that differ in ids are not identical."""
    reaction = model.reactions.get_by_id("PGI")
    old_reaction = reaction.copy()
    reaction.id = "PGI2"
    equivalent, comparison = compare_reaction_state(reaction, old_reaction)
    assert not equivalent
    assert comparison["same"] == set(reaction.__getstate__().keys()).difference({"_id"})
    assert comparison["modified"] == {"_id": ("PGI2", "PGI")}
    assert comparison["added"] == set()
    assert comparison["removed"] == set()


def test_reaction_different_names(model: Model) -> None:
    """Test that reactions that differ in names are not identical."""
    reaction = model.reactions.get_by_id("PGI")
    old_reaction = reaction.copy()
    reaction.name = reaction.name + " "
    equivalent, comparison = compare_reaction_state(reaction, old_reaction)
    assert not equivalent
    assert comparison["same"] == set(reaction.__getstate__().keys()).difference(
        {"name"}
    )
    assert comparison["modified"] == {
        "name": (old_reaction.name + " ", old_reaction.name)
    }
    assert comparison["added"] == set()
    assert comparison["removed"] == set()
    reaction.name = None
    equivalent, comparison = compare_reaction_state(reaction, old_reaction)
    assert not equivalent
    assert comparison["same"] == set(reaction.__getstate__().keys()).difference(
        {"name"}
    )
    assert comparison["modified"] == {"name": (None, old_reaction.name)}
    assert comparison["added"] == set()
    assert comparison["removed"] == set()
    reaction.name = ""
    equivalent, comparison = compare_reaction_state(reaction, old_reaction)
    assert not equivalent
    assert comparison["same"] == set(reaction.__getstate__().keys()).difference(
        {"name"}
    )
    assert comparison["modified"] == {"name": ("", old_reaction.name)}
    assert comparison["added"] == set()
    assert comparison["removed"] == set()
    reaction.name = "Test"
    equivalent, comparison = compare_reaction_state(reaction, old_reaction)
    assert not equivalent
    assert comparison["same"] == set(reaction.__getstate__().keys()).difference(
        {"name"}
    )
    assert comparison["modified"] == {"name": ("Test", old_reaction.name)}
    assert comparison["added"] == set()
    assert comparison["removed"] == set()


def test_reaction_gpr_modification(model: Model) -> None:
    """Test reactions are not equivalent after GPR/gene rule manipulations."""
    reaction = model.reactions.get_by_id("PGI")
    old_reaction = reaction.copy()
    old_gene = list(reaction.genes)[0]

    # Add an existing 'gene' to the GPR
    reaction.gene_reaction_rule = "s0001"
    equivalent, comparison = compare_reaction_state(reaction, old_reaction)
    assert not equivalent
    assert comparison["modified"]["_genes"] == ({"s0001"}, {old_gene.id})
    assert comparison["modified"]["_gpr"] == (reaction.gpr, old_reaction.gpr)
    assert comparison["same"] == set(reaction.__getstate__().keys()).difference(
        {"_gpr", "_genes"}
    )

    old_reaction = reaction.copy()
    old_reaction.gpr = GPR()
    equivalent, comparison = compare_reaction_state(reaction, old_reaction)
    assert not equivalent
    assert comparison["modified"]["_genes"] == ({"s0001"}, set())
    assert comparison["modified"]["_gpr"] == (reaction.gpr, GPR())
    assert comparison["same"] == set(reaction.__getstate__().keys()).difference(
        {"_gpr", "_genes"}
    )


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
    assert comparison["same"] == set(acald_reaction.__getstate__().keys()).difference(
        {"_lower_bound", "_upper_bound"}
    )
    assert comparison["modified"] == {
        "_lower_bound": (-1000, -1100),
        "_upper_bound": (1000, -1100),
    }
    assert comparison["added"] == set()
    assert comparison["removed"] == set()
    new_reaction = acald_reaction.copy()
    new_reaction.bounds = (
        acald_reaction.upper_bound + 100,
        acald_reaction.upper_bound + 100,
    )
    equivalent, comparison = compare_reaction_state(acald_reaction, new_reaction)
    assert not equivalent
    assert comparison["same"] == set(acald_reaction.__getstate__().keys()).difference(
        {"_lower_bound", "_upper_bound"}
    )
    assert comparison["modified"] == {
        "_lower_bound": (-1000, 1100),
        "_upper_bound": (1000, 1100),
    }
    assert comparison["added"] == set()
    assert comparison["removed"] == set()
    new_reaction = acald_reaction.copy()
    new_reaction.lower_bound = -100
    equivalent, comparison = compare_reaction_state(acald_reaction, new_reaction)
    assert not equivalent
    assert comparison["same"] == set(acald_reaction.__getstate__().keys()).difference(
        {"_lower_bound"}
    )
    assert comparison["modified"] == {"_lower_bound": (-1000, -100)}
    assert comparison["added"] == set()
    assert comparison["removed"] == set()
    new_reaction = acald_reaction.copy()
    new_reaction.upper_bound = 100
    equivalent, comparison = compare_reaction_state(acald_reaction, new_reaction)
    assert not equivalent
    assert comparison["same"] == set(acald_reaction.__getstate__().keys()).difference(
        {"_upper_bound"}
    )
    assert comparison["modified"] == {"_upper_bound": (1000, 100)}
    assert comparison["added"] == set()
    assert comparison["removed"] == set()
    new_reaction = acald_reaction.copy()
    new_reaction.knock_out()
    equivalent, comparison = compare_reaction_state(acald_reaction, new_reaction)
    assert not equivalent
    assert comparison["same"] == set(acald_reaction.__getstate__().keys()).difference(
        {"_lower_bound", "_upper_bound"}
    )
    assert comparison["modified"] == {
        "_lower_bound": (-1000, 0),
        "_upper_bound": (1000, 0),
    }
    assert comparison["added"] == set()
    assert comparison["removed"] == set()
    new_reaction = acald_reaction.copy()
    new_reaction.bounds = (0, 0)
    equivalent, comparison = compare_reaction_state(acald_reaction, new_reaction)
    assert not equivalent
    assert comparison["same"] == set(acald_reaction.__getstate__().keys()).difference(
        {"_lower_bound", "_upper_bound"}
    )
    assert comparison["modified"] == {
        "_lower_bound": (-1000, 0),
        "_upper_bound": (1000, 0),
    }
    assert comparison["added"] == set()
    assert comparison["removed"] == set()


def test_compare_reaction_different_subsystem(model: Model) -> None:
    """Test that differences in subsystem are picked up by compare function."""
    acald_reaction = model.reactions.ACALD
    new_reaction = acald_reaction.copy()
    new_reaction.subsystem = "Test"
    equivalent, comparison = compare_reaction_state(acald_reaction, new_reaction)
    assert not equivalent
    assert comparison["same"] == set(acald_reaction.__getstate__().keys()).difference(
        {"subsystem"}
    )
    assert comparison["modified"] == {"subsystem": ("", "Test")}
    new_reaction = acald_reaction.copy()
    new_reaction.subsystem = None
    equivalent, comparison = compare_reaction_state(acald_reaction, new_reaction)
    assert not equivalent
    assert comparison["same"] == set(acald_reaction.__getstate__().keys()).difference(
        {"subsystem"}
    )
    assert comparison["modified"] == {"subsystem": ("", None)}


def test_add_metabolite_comparison(model: Model) -> None:
    """Test metabolite addition to a reaction is not equivalent to original reaction."""
    with model:
        with model:
            reaction = model.reactions.get_by_id("PGI")
            old_reaction = reaction.copy()
            reaction.add_metabolites({model.metabolites[0]: 1})
            equivalent, comparison = compare_reaction_state(old_reaction, reaction)
            assert not equivalent
            assert comparison["same"] == set(reaction.__getstate__().keys()).difference(
                {"_metabolites"}
            )
            old_metabolites = {k.id: v for k, v in old_reaction.metabolites.items()}
            reaction_metabolites = old_metabolites.copy()
            reaction_metabolites.update({model.metabolites[0].id: 1})
            assert comparison["modified"] == {
                "_metabolites": (old_metabolites, reaction_metabolites)
            }

            fake_metabolite = Metabolite("fake")
            reaction.add_metabolites({fake_metabolite: 1})
            equivalent, comparison = compare_reaction_state(old_reaction, reaction)
            assert not equivalent
            assert comparison["same"] == set(reaction.__getstate__().keys()).difference(
                {"_metabolites"}
            )
            reaction_metabolites.update({"fake": 1})
            assert comparison["modified"] == {
                "_metabolites": (old_metabolites, reaction_metabolites)
            }

    equivalent, comparison = compare_reaction_state(old_reaction, reaction)
    assert equivalent
    assert comparison["same"] == set(reaction.__getstate__().keys())

    # Test adding by string
    with model:
        reaction.add_metabolites({"g6p_c": -1})  # already in reaction
        reaction_metabolites = {m.id: v for m, v in old_reaction.metabolites.items()}
        reaction_metabolites["g6p_c"] = -2
        equivalent, comparison = compare_reaction_state(old_reaction, reaction)
        assert not equivalent
        assert comparison["modified"] == {
            "_metabolites": (old_metabolites, reaction_metabolites)
        }
        reaction.add_metabolites({"h_c": 1})
        reaction_metabolites["h_c"] = 1
        equivalent, comparison = compare_reaction_state(old_reaction, reaction)
        assert not equivalent
        assert comparison["modified"] == {
            "_metabolites": (old_metabolites, reaction_metabolites)
        }

    equivalent, _ = compare_reaction_state(old_reaction, reaction)
    assert equivalent

    # Test combine=False
    reaction = model.reactions.get_by_id("ATPM")
    old_reaction = reaction.copy()
    old_metabolites = {k.id: v for k, v in old_reaction.metabolites.items()}
    reaction_metabolites = old_metabolites.copy()
    reaction_metabolites["h2o_c"] = 2.5
    with model:
        reaction.add_metabolites({"h2o_c": 2.5}, combine=False)
        equivalent, comparison = compare_reaction_state(old_reaction, reaction)
        assert not equivalent
        assert comparison["modified"] == {
            "_metabolites": (old_metabolites, reaction_metabolites)
        }

    # Test adding to a new Reaction
    reaction = Reaction("test")
    old_reaction = reaction.copy()
    reaction.add_metabolites({Metabolite("test_met"): -1})
    equivalent, comparison = compare_reaction_state(old_reaction, reaction)
    assert not equivalent
    assert comparison["modified"] == {"_metabolites": ({}, {"test_met": -1})}


def test_iadd_reaction_comparison(model: Model) -> None:
    """Test in-place addition of reaction is correctly identified by compare."""
    PGI = model.reactions.PGI
    PGI_copy = PGI.copy()
    EX_h2o = model.reactions.EX_h2o_e
    PGI += EX_h2o
    PGI_metabolites = {k.id: v for k, v in PGI_copy.metabolites.items()}
    PGI_copy_metabolites = PGI_metabolites.copy()
    PGI_copy_metabolites[model.metabolites.h2o_e.id] = -1.0

    equivalent, comparison = compare_reaction_state(PGI_copy, PGI)
    assert not equivalent
    assert comparison["same"] == set(PGI.__getstate__().keys()).difference(
        {"_metabolites"}
    )
    assert comparison["modified"] == {
        "_metabolites": (PGI_metabolites, PGI_copy_metabolites)
    }
    # Add a reaction not in the model
    new_reaction = Reaction("test")
    new_reaction.add_metabolites({Metabolite("A"): -1, Metabolite("B"): 1})
    PGI += new_reaction
    PGI_copy_metabolites["A"] = -1
    PGI_copy_metabolites["B"] = 1
    equivalent, comparison = compare_reaction_state(PGI_copy, PGI)
    assert not equivalent
    assert comparison["same"] == set(PGI.__getstate__().keys()).difference(
        {"_metabolites"}
    )
    assert comparison["modified"] == {
        "_metabolites": (PGI_metabolites, PGI_copy_metabolites)
    }
    # Combine two GPRs
    ACKr_copy = model.reactions.ACKr.copy()
    ACKr = model.reactions.ACKr
    model.reactions.ACKr += model.reactions.ACONTa
    expected_genes = {g.id for g in ACKr.genes}
    expected_rule = "(b2296 or b3115 or b1849) and (b0118 or b1276)"
    ACKr_metabolites = {m.id: stoic for m, stoic in ACKr_copy.metabolites.items()}
    expected_metabolites = ACKr_metabolites.copy()
    expected_metabolites.update(
        {m.id: stoic for m, stoic in model.reactions.ACONTa.metabolites.items()}
    )
    equivalent, comparison = compare_reaction_state(ACKr_copy, ACKr)
    assert not equivalent
    assert comparison["same"] == set(PGI.__getstate__().keys()).difference(
        {"_metabolites", "_gpr", "_genes"}
    )
    assert comparison["modified"]["_genes"] == (
        {g.id for g in ACKr_copy.genes},
        expected_genes,
    )
    assert comparison["modified"]["_gpr"] == (
        ACKr_copy.gpr,
        GPR().from_string(expected_rule),
    )
    assert comparison["modified"]["_metabolites"] == (
        ACKr_metabolites,
        expected_metabolites,
    )
    assert comparison["modified"] == {
        "_metabolites": (ACKr_metabolites, expected_metabolites),
        "_genes": ({g.id for g in ACKr_copy.genes}, expected_genes),
        "_gpr": (ACKr_copy.gpr, GPR().from_string(expected_rule)),
    }


def test_mul_reaction_comparison(model: Model) -> None:
    """Test scalar multiplication of factors with a reaction."""
    new = model.reactions.PGI * 2
    PGI = model.reactions.PGI
    equivalent, comparison = compare_reaction_state(PGI, new)
    assert not equivalent
    assert comparison["same"] == set(PGI.__getstate__().keys()).difference(
        {"_metabolites"}
    )
    PGI_metabolites = {m.id: stoic for m, stoic in PGI.metabolites.items()}
    new_metabolites = {m.id: stoic * 2 for m, stoic in PGI.metabolites.items()}
    assert comparison["modified"] == {
        "_metabolites": (PGI_metabolites, new_metabolites)
    }


def test_sub_reaction_comparison(model: Model) -> None:
    """Test reaction subtraction is picked up by comparison function."""
    new = model.reactions.PGI - model.reactions.EX_h2o_e
    PGI = model.reactions.PGI
    PGI_metabolites = {m.id: stoic for m, stoic in PGI.metabolites.items()}
    new_metabolites = PGI_metabolites.copy()
    new_metabolites[model.metabolites.h2o_e.id] = 1.0
    equivalent, comparison = compare_reaction_state(PGI, new)
    assert not equivalent
    assert comparison["same"] == set(PGI.__getstate__().keys()).difference(
        {"_metabolites"}
    )
    assert not equivalent
    assert comparison["modified"] == {
        "_metabolites": (PGI_metabolites, new_metabolites)
    }


def test_add_metabolites_combine_true_reaction_comparison(model: Model) -> None:
    """Test metabolite addition to reaction (combine = True) with comparison."""
    test_metabolite = Metabolite("test")
    for reaction in model.reactions:
        old_reaction = reaction.copy()
        reaction.add_metabolites({test_metabolite: -66}, combine=True)
        reaction_metabolites = {
            m.id: stoic for m, stoic in reaction.metabolites.items()
        }
        old_reaction_metabolites = reaction_metabolites.copy()
        old_reaction_metabolites.pop("test")
        equivalent, comparison = compare_reaction_state(old_reaction, reaction)
        assert not equivalent
        assert comparison["same"] == set(reaction.__getstate__().keys()).difference(
            {"_metabolites"}
        )
        assert comparison["modified"] == {
            "_metabolites": (old_reaction_metabolites, reaction_metabolites)
        }
        already_included_metabolite = list(reaction.metabolites.keys())[0]
        previous_coefficient = reaction.get_coefficient(already_included_metabolite.id)
        old_reaction = reaction.copy()
        reaction.add_metabolites({already_included_metabolite: 10}, combine=True)
        reaction_metabolites = {
            m.id: stoic for m, stoic in reaction.metabolites.items()
        }
        old_reaction_metabolites = reaction_metabolites.copy()
        old_reaction_metabolites[already_included_metabolite.id] = previous_coefficient
        equivalent, comparison = compare_reaction_state(old_reaction, reaction)
        assert comparison["same"] == set(reaction.__getstate__().keys()).difference(
            {"_metabolites"}
        )
        assert comparison["modified"] == {
            "_metabolites": (old_reaction_metabolites, reaction_metabolites)
        }


def test_reaction_annotation_comparison(model: Model) -> None:
    """Test that changes in annotation are picked up by comparison.

    Parameters
    ----------
    model: cobra.Model


    """
    PGI = model.reactions.PGI
    PGI_copy = PGI.copy()
    PGI_copy_annotation = PGI_copy.annotation
    with model:
        PGI.annotation = {}
        equivalent, comparison = compare_reaction_state(PGI_copy, PGI)
        assert not equivalent
        assert comparison["same"] == set(PGI.__getstate__().keys()).difference(
            {"_annotation"}
        )
        assert comparison["modified"] == {"_annotation": (PGI_copy_annotation, {})}
    with model:
        PGI.annotation = {k: v + " " for k, v in PGI.annotation.items()}
        equivalent, comparison = compare_reaction_state(PGI_copy, PGI)
        assert not equivalent
        assert comparison["same"] == set(PGI.__getstate__().keys()).difference(
            {"_annotation"}
        )
        assert comparison["modified"] == {
            "_annotation": (
                PGI_copy_annotation,
                {k: v + " " for k, v in PGI_copy.annotation.items()},
            )
        }
    with model:
        PGI.annotation["bigg.reaction"] = "Test"
        equivalent, comparison = compare_reaction_state(PGI_copy, PGI)
        PGI_annotation = PGI_copy_annotation
        PGI_annotation["bigg.reaction"] = "Test"
        assert not equivalent
        assert comparison["same"] == set(PGI.__getstate__().keys()).difference(
            {"_annotation"}
        )
        assert comparison["modified"] == {
            "_annotation": (PGI_copy_annotation, PGI_annotation)
        }


def test_reaction_notes_comparison(model: Model) -> None:
    """Test that notes in reaction can be picked up by comparison function."""
    PGI = model.reactions.PGI
    PGI_copy = PGI.copy()
    PGI.notes = {"Note": "Test"}
    equivalent, comparison = compare_reaction_state(PGI_copy, PGI)
    assert not equivalent
    assert comparison["same"] == set(PGI.__getstate__().keys()).difference({"notes"})
    assert comparison["modified"] == {"notes": ({}, {"Note": "Test"})}
    PGI_copy = PGI.copy()
    PGI.notes = {"Note": "Test "}
    equivalent, comparison = compare_reaction_state(PGI_copy, PGI)
    assert not equivalent
    assert comparison["same"] == set(PGI.__getstate__().keys()).difference({"notes"})
    assert comparison["modified"] == {"notes": ({"Note": "Test"}, {"Note": "Test "})}
    PGI.notes = {"Note": "test"}
    equivalent, comparison = compare_reaction_state(PGI_copy, PGI)
    assert not equivalent
    assert comparison["same"] == set(PGI.__getstate__().keys()).difference({"notes"})
    assert comparison["modified"] == {"notes": ({"Note": "Test"}, {"Note": "test"})}
    PGI.notes = {"note": "Test"}
    equivalent, comparison = compare_reaction_state(PGI_copy, PGI)
    assert not equivalent
    assert comparison["same"] == set(PGI.__getstate__().keys()).difference({"notes"})
    assert comparison["modified"] == {"notes": ({"Note": "Test"}, {"note": "Test"})}
    PGI.notes = {"Note": "Test", "secondNote": "test"}
    equivalent, comparison = compare_reaction_state(PGI_copy, PGI)
    assert not equivalent
    assert comparison["same"] == set(PGI.__getstate__().keys()).difference({"notes"})
    assert comparison["modified"] == {
        "notes": (
            {"Note": "Test"},
            {"Note": "Test", "secondNote": "test"},
        )
    }


def test_reaction_comparison_ignore_keys(model: Model) -> None:
    """Test that the ignore_keys field in reaction comparison works as expected."""
    PGI = model.reactions.get_by_id("PGI")
    PGI_copy = PGI.copy()
    PGI.blah = None
    equivalent, comparison = compare_reaction_state(PGI, PGI_copy, ignore_keys={"blah"})
    assert equivalent
    assert comparison["same"] == PGI_copy.__getstate__().keys()
    assert comparison["modified"] == {}
    assert comparison["added"] == set()
    assert comparison["removed"] == set()
    assert not hasattr(PGI_copy, "blah")
    assert hasattr(PGI, "blah")
    PGI.__dict__.pop("blah")

    PGI.id = "PGI2"
    equivalent, comparison = compare_reaction_state(PGI, PGI_copy, ignore_keys={"_id"})
    assert equivalent
    assert comparison["same"] == set(PGI_copy.__getstate__().keys()).difference({"_id"})
    assert comparison["modified"] == {}
    assert comparison["added"] == set()
    assert comparison["removed"] == set()
    assert PGI.id != PGI_copy.id
    PGI.id = PGI_copy.id

    PGI.name = PGI.name + " "
    equivalent, comparison = compare_reaction_state(PGI, PGI_copy, ignore_keys={"name"})
    assert equivalent
    assert comparison["same"] == set(PGI_copy.__getstate__().keys()).difference(
        {"name"}
    )
    assert comparison["modified"] == {}
    assert comparison["added"] == set()
    assert comparison["removed"] == set()
    assert PGI.name != PGI_copy.name
    PGI.name = PGI_copy.name

    with model:
        PGI.gene_reaction_rule = "s0001"
        equivalent, comparison = compare_reaction_state(
            PGI, PGI_copy, ignore_keys={"_genes", "_gpr"}
        )
        assert equivalent
        assert comparison["same"] == set(PGI_copy.__getstate__().keys()).difference(
            {"_genes", "_gpr"}
        )
        assert comparison["modified"] == {}
        assert comparison["added"] == set()
        assert comparison["removed"] == set()
        assert PGI.gpr != PGI_copy.gpr

    with model:
        PGI.bounds = (
            PGI_copy.lower_bound - 100,
            PGI_copy.lower_bound - 100,
        )
        equivalent, comparison = compare_reaction_state(
            PGI, PGI_copy, ignore_keys={"_lower_bound", "_upper_bound"}
        )
        assert equivalent
        assert comparison["same"] == set(PGI_copy.__getstate__().keys()).difference(
            {"_lower_bound", "_upper_bound"}
        )
        assert comparison["modified"] == {}
        assert comparison["added"] == set()
        assert comparison["removed"] == set()
        assert PGI.bounds != PGI_copy.bounds

    PGI.subsystem = "Test"
    equivalent, comparison = compare_reaction_state(
        PGI, PGI_copy, ignore_keys={"subsystem"}
    )
    assert equivalent
    assert comparison["same"] == set(PGI_copy.__getstate__().keys()).difference(
        {"subsystem"}
    )
    assert comparison["modified"] == {}
    assert comparison["added"] == set()
    assert comparison["removed"] == set()
    assert PGI.subsystem != PGI_copy.subsystem
    PGI.subsystem = PGI_copy.subsystem

    with model:
        PGI.add_metabolites({model.metabolites[0]: 1})
        equivalent, comparison = compare_reaction_state(
            PGI, PGI_copy, ignore_keys={"_metabolites"}
        )
        assert equivalent
        assert comparison["same"] == set(PGI_copy.__getstate__().keys()).difference(
            {"_metabolites"}
        )
        assert comparison["modified"] == {}
        assert comparison["added"] == set()
        assert comparison["removed"] == set()
        assert PGI.metabolites != PGI_copy.metabolites
    with model:
        PGI.annotation = {}
        equivalent, comparison = compare_reaction_state(
            PGI, PGI_copy, ignore_keys={"_annotation"}
        )
        assert equivalent
        assert comparison["same"] == set(PGI_copy.__getstate__().keys()).difference(
            {"_annotation"}
        )
        assert comparison["modified"] == {}
        assert comparison["added"] == set()
        assert comparison["removed"] == set()
        assert PGI.annotation != PGI_copy.annotation
    with model:
        PGI.notes = {"Note": "Test"}
        equivalent, comparison = compare_reaction_state(
            PGI, PGI_copy, ignore_keys={"notes"}
        )
        assert equivalent
        assert comparison["same"] == set(PGI_copy.__getstate__().keys()).difference(
            {"notes"}
        )
        assert comparison["modified"] == {}
        assert comparison["added"] == set()
        assert comparison["removed"] == set()
        assert PGI.notes != PGI_copy.notes


def test_metabolite_copies_are_equivalent(model: Model) -> None:
    """Test that a copy of a metabolite will return True with compare functions."""
    NADH = model.metabolites.get_by_id("nadh_c")
    NADH_copy = NADH.copy()
    equivalent, comparison = compare_state(NADH, NADH_copy)
    assert equivalent
    assert comparison["same"] == NADH.__getstate__().keys()
    assert comparison["modified"] == {}
    assert comparison["added"] == set()
    assert comparison["removed"] == set()


def test_metabolite_with_added_field_is_different(model: Model) -> None:
    """Test that metabolites with added fields are picked up by comparison."""
    NADH = model.metabolites.get_by_id("nadh_c")
    NADH_copy = NADH.copy()
    NADH.blah = None
    equivalent, comparison = compare_state(NADH, NADH_copy)
    assert not equivalent
    assert comparison["same"] == NADH_copy.__getstate__().keys()
    assert comparison["modified"] == {}
    assert comparison["added"] == {"blah"}
    assert comparison["removed"] == set()
    equivalent, comparison = compare_state(NADH_copy, NADH)
    assert not equivalent
    assert comparison["same"] == NADH_copy.__getstate__().keys()
    assert comparison["modified"] == {}
    assert comparison["added"] == set()
    assert comparison["removed"] == {"blah"}


@pytest.mark.parametrize("field_name", ["_id", "name", "formula", "compartment"])
def test_metabolite_comparison_different_string_fields(
    model: Model, field_name: str
) -> None:
    """Test that metabolites that differ in string fields are not identical.

    This function will test id (_id), name, formula, compartment.

    Parameters
    ----------
    model: cobra.Model
        Model to take metabolites from
    field_name: str
        Which field to test.

    """
    NADH = model.metabolites.get_by_id("nadh_c")
    NADH_copy = NADH.copy()
    NADH.__setattr__(field_name, NADH.__getattribute__(field_name) + " ")
    equivalent, comparison = compare_state(NADH, NADH_copy)
    assert not equivalent
    assert comparison["same"] == set(NADH.__getstate__().keys()).difference(
        {field_name}
    )
    assert comparison["modified"] == {
        field_name: (
            NADH_copy.__getattribute__(field_name) + " ",
            NADH_copy.__getattribute__(field_name),
        )
    }
    assert comparison["added"] == set()
    assert comparison["removed"] == set()
    NADH.__setattr__(field_name, None)
    equivalent, comparison = compare_state(NADH, NADH_copy)
    assert not equivalent
    assert comparison["same"] == set(NADH.__getstate__().keys()).difference(
        {field_name}
    )
    assert comparison["modified"] == {
        field_name: (None, NADH_copy.__getattribute__(field_name))
    }
    assert comparison["added"] == set()
    assert comparison["removed"] == set()
    NADH.__setattr__(field_name, "")
    equivalent, comparison = compare_state(NADH, NADH_copy)
    assert not equivalent
    assert comparison["same"] == set(NADH.__getstate__().keys()).difference(
        {field_name}
    )
    assert comparison["modified"] == {
        field_name: ("", NADH_copy.__getattribute__(field_name))
    }
    assert comparison["added"] == set()
    assert comparison["removed"] == set()
    NADH.__setattr__(field_name, "Test")
    equivalent, comparison = compare_state(NADH, NADH_copy)
    assert not equivalent
    assert comparison["same"] == set(NADH.__getstate__().keys()).difference(
        {field_name}
    )
    assert comparison["modified"] == {
        field_name: ("Test", NADH_copy.__getattribute__(field_name))
    }
    assert comparison["added"] == set()
    assert comparison["removed"] == set()
    NADH.__setattr__(field_name, "C21H27N7O14P")
    equivalent, comparison = compare_state(NADH, NADH_copy)
    assert not equivalent
    assert comparison["same"] == set(NADH.__getstate__().keys()).difference(
        {field_name}
    )
    assert comparison["modified"] == {
        field_name: ("C21H27N7O14P", NADH_copy.__getattribute__(field_name))
    }
    assert comparison["added"] == set()
    assert comparison["removed"] == set()
    NADH.__setattr__(field_name, "e")
    equivalent, comparison = compare_state(NADH, NADH_copy)
    assert not equivalent
    assert comparison["same"] == set(NADH.__getstate__().keys()).difference(
        {field_name}
    )
    assert comparison["modified"] == {
        field_name: ("e", NADH_copy.__getattribute__(field_name))
    }
    assert comparison["added"] == set()
    assert comparison["removed"] == set()


def test_metabolite_charge_comparison(model: Model) -> None:
    """Test that metabolites with different charge are picked up by comparison."""
    NADH = model.metabolites.get_by_id("nadh_c")
    NADH_copy = NADH.copy()
    NADH.charge = 0
    equivalent, comparison = compare_state(NADH, NADH_copy)
    assert not equivalent
    assert comparison["same"] == set(NADH.__getstate__().keys()).difference({"charge"})
    assert comparison["modified"] == {"charge": (0, NADH_copy.charge)}
    assert comparison["added"] == set()
    assert comparison["removed"] == set()
    NADH.charge = None
    equivalent, comparison = compare_state(NADH, NADH_copy)
    assert not equivalent
    assert comparison["same"] == set(NADH.__getstate__().keys()).difference({"charge"})
    assert comparison["modified"] == {"charge": (None, NADH_copy.charge)}
    assert comparison["added"] == set()
    assert comparison["removed"] == set()
    NADH.charge = 2
    equivalent, comparison = compare_state(NADH, NADH_copy)
    assert not equivalent
    assert comparison["same"] == set(NADH.__getstate__().keys()).difference({"charge"})
    assert comparison["modified"] == {"charge": (2, NADH_copy.charge)}
    assert comparison["added"] == set()
    assert comparison["removed"] == set()


def test_metabolite_annotation_comparison(model: Model) -> None:
    """Test that changes in metabolite annotation are picked up by comparison."""
    NADH = model.metabolites.get_by_id("nadh_c")
    NADH_copy = NADH.copy()
    NADH_copy_annotation = NADH_copy.annotation
    with model:
        NADH.annotation = {}
        equivalent, comparison = compare_state(NADH_copy, NADH)
        assert not equivalent
        assert comparison["same"] == set(NADH.__getstate__().keys()).difference(
            {"_annotation"}
        )
        assert comparison["modified"] == {"_annotation": (NADH_copy_annotation, {})}
        assert comparison["added"] == set()
        assert comparison["removed"] == set()
    with model:
        NADH_annotation = {k: v for k, v in NADH_copy_annotation.items()}
        NADH_annotation["bigg.metabolite"] = "Test"
        NADH.annotation = NADH_annotation
        equivalent, comparison = compare_state(NADH_copy, NADH)
        assert not equivalent
        assert comparison["same"] == set(NADH.__getstate__().keys()).difference(
            {"_annotation"}
        )
        assert comparison["modified"] == {
            "_annotation": (NADH_copy_annotation, NADH_annotation)
        }
        assert comparison["added"] == set()
        assert comparison["removed"] == set()
    with model:
        NADH_annotation = {k: v for k, v in NADH_copy_annotation.items()}
        NADH_annotation.pop("biocyc")
        NADH.annotation = NADH_annotation
        equivalent, comparison = compare_state(NADH_copy, NADH)
        assert not equivalent
        assert comparison["same"] == set(NADH.__getstate__().keys()).difference(
            {"_annotation"}
        )
        assert comparison["modified"] == {
            "_annotation": (NADH_copy_annotation, NADH_annotation)
        }
        assert comparison["added"] == set()
        assert comparison["removed"] == set()
    with model:
        NADH_annotation = {k: v for k, v in NADH_copy_annotation.items()}
        NADH_annotation["test"] = "test"
        NADH.annotation = NADH_annotation
        equivalent, comparison = compare_state(NADH_copy, NADH)
        assert not equivalent
        assert comparison["same"] == set(NADH.__getstate__().keys()).difference(
            {"_annotation"}
        )
        assert comparison["modified"] == {
            "_annotation": (NADH_copy_annotation, NADH_annotation)
        }
        assert comparison["added"] == set()
        assert comparison["removed"] == set()
    with model:
        NADH_annotation = {
            "bigg.metabolite": "nadh ",
            "biocyc": "NADH ",
            "cas": ["58-68-4 "],
            "chebi": [
                "CHEBI:13395 ",
                "CHEBI:21902 ",
                "CHEBI:16908 ",
                "CHEBI:7423 ",
                "CHEBI:44216 ",
                "CHEBI:57945 ",
                "CHEBI:13396 ",
            ],
            "hmdb": "HMDB01487 ",
            "kegg.compound": "C00004 ",
            "pubchem.substance": "3306 ",
            "reactome": [
                "REACT_192305 ",
                "REACT_73473 ",
                "REACT_194697 ",
                "REACT_29362 ",
            ],
            "seed.compound": "cpd00004 ",
            "unipathway.compound": "UPC00004 ",
        }
        NADH.annotation = NADH_annotation
        equivalent, comparison = compare_state(NADH_copy, NADH)
        assert not equivalent
        assert comparison["same"] == set(NADH.__getstate__().keys()).difference(
            {"_annotation"}
        )
        assert comparison["modified"] == {
            "_annotation": (NADH_copy_annotation, NADH_annotation)
        }
        assert comparison["added"] == set()
        assert comparison["removed"] == set()


def test_metabolite_notes_comparison(model: Model) -> None:
    """Test that notes in Metabolite can be picked up by comparison function."""
    NADH = model.metabolites.get_by_id("nadh_c")
    NADH_copy = NADH.copy()
    NADH.notes = {"Note": "Test"}
    equivalent, comparison = compare_state(NADH_copy, NADH)
    assert not equivalent
    assert comparison["same"] == set(NADH.__getstate__().keys()).difference({"notes"})
    assert comparison["modified"] == {"notes": ({}, {"Note": "Test"})}
    NADH_copy = NADH.copy()
    NADH.notes = {"Note": "Test "}
    equivalent, comparison = compare_state(NADH_copy, NADH)
    assert not equivalent
    assert comparison["same"] == set(NADH.__getstate__().keys()).difference({"notes"})
    assert comparison["modified"] == {"notes": ({"Note": "Test"}, {"Note": "Test "})}
    NADH.notes = {"Note": "test"}
    equivalent, comparison = compare_state(NADH_copy, NADH)
    assert not equivalent
    assert comparison["same"] == set(NADH.__getstate__().keys()).difference({"notes"})
    assert comparison["modified"] == {"notes": ({"Note": "Test"}, {"Note": "test"})}
    NADH.notes = {"note": "Test"}
    equivalent, comparison = compare_state(NADH_copy, NADH)
    assert not equivalent
    assert comparison["same"] == set(NADH.__getstate__().keys()).difference({"notes"})
    assert comparison["modified"] == {"notes": ({"Note": "Test"}, {"note": "Test"})}
    NADH.notes = {"Note": "Test", "secondNote": "test"}
    equivalent, comparison = compare_state(NADH_copy, NADH)
    assert not equivalent
    assert comparison["same"] == set(NADH.__getstate__().keys()).difference({"notes"})
    assert comparison["modified"] == {
        "notes": (
            {"Note": "Test"},
            {"Note": "Test", "secondNote": "test"},
        )
    }


def test_metabolite_comparison_ignore_keys(model: Model) -> None:
    """Test that the ignore_keys field in Metabolite comparison works as expected."""
    NADH = model.metabolites.get_by_id("nadh_c")
    NADH_copy = NADH.copy()
    NADH.blah = None
    equivalent, comparison = compare_state(NADH, NADH_copy, ignore_keys={"blah"})
    assert equivalent
    assert comparison["same"] == NADH_copy.__getstate__().keys()
    assert comparison["modified"] == {}
    assert comparison["added"] == set()
    assert comparison["removed"] == set()
    NADH.__dict__.pop("blah")

    NADH.charge = 0
    equivalent, comparison = compare_state(NADH, NADH_copy, {"charge"})
    assert equivalent
    assert comparison["same"] == set(NADH.__getstate__().keys()).difference({"charge"})
    assert comparison["modified"] == {}
    assert comparison["added"] == set()
    assert comparison["removed"] == set()
    assert NADH.charge != NADH_copy.charge
    NADH.charge = NADH_copy.charge

    NADH_copy_annotation = NADH_copy.annotation
    with model:
        NADH_annotation = {k: v for k, v in NADH_copy_annotation.items()}
        NADH_annotation["bigg.metabolite"] = "Test"
        NADH.annotation = NADH_annotation
        equivalent, comparison = compare_state(
            NADH_copy, NADH, ignore_keys={"_annotation"}
        )
        assert equivalent
        assert comparison["same"] == set(NADH.__getstate__().keys()).difference(
            {"_annotation"}
        )
        assert comparison["modified"] == {}
        assert comparison["added"] == set()
        assert comparison["removed"] == set()
        assert NADH.annotation != NADH_copy.annotation
    with model:
        NADH_annotation = {k: v for k, v in NADH_copy_annotation.items()}
        NADH_annotation.pop("biocyc")
        NADH.annotation = NADH_annotation
        equivalent, comparison = compare_state(
            NADH_copy, NADH, ignore_keys={"_annotation"}
        )
        assert equivalent
        assert comparison["same"] == set(NADH.__getstate__().keys()).difference(
            {"_annotation"}
        )
        assert comparison["modified"] == {}
        assert comparison["added"] == set()
        assert comparison["removed"] == set()
        assert NADH.annotation != NADH_copy.annotation

    NADH.notes = {"Note": "Test"}
    equivalent, comparison = compare_state(NADH_copy, NADH, ignore_keys={"notes"})
    assert equivalent
    assert comparison["same"] == set(NADH.__getstate__().keys()).difference({"notes"})
    assert comparison["modified"] == {}
    assert comparison["added"] == set()
    assert comparison["removed"] == set()
    assert NADH.notes != NADH_copy.notes


@pytest.mark.parametrize("field_name", ["_id", "name", "formula", "compartment"])
def test_metabolite_comparison_ignore_keys_different_string_fields(
    model: Model, field_name: str
) -> None:
    """Test that ignore keys works on string fields in metaoblites.

    This function will test id (_id), name, formula, compartment.

    Parameters
    ----------
    model: cobra.Model
        Model to take metabolites from
    field_name: str
        Which field to test.

    """
    NADH = model.metabolites.get_by_id("nadh_c")
    NADH_copy = NADH.copy()
    NADH.__setattr__(field_name, NADH.__getattribute__(field_name) + " ")
    equivalent, comparison = compare_state(NADH, NADH_copy, ignore_keys={field_name})
    assert equivalent
    assert comparison["same"] == set(NADH.__getstate__().keys()).difference(
        {field_name}
    )
    assert comparison["modified"] == {}
    assert comparison["added"] == set()
    assert comparison["removed"] == set()
    assert NADH.__getattribute__(field_name) != NADH_copy.__getattribute__(field_name)
    NADH.__setattr__(field_name, None)
    equivalent, comparison = compare_state(NADH, NADH_copy, ignore_keys={field_name})
    assert equivalent
    assert comparison["same"] == set(NADH.__getstate__().keys()).difference(
        {field_name}
    )
    assert comparison["modified"] == {}
    assert comparison["added"] == set()
    assert comparison["removed"] == set()
    assert NADH.__getattribute__(field_name) != NADH_copy.__getattribute__(field_name)
    NADH.__setattr__(field_name, "")
    equivalent, comparison = compare_state(NADH, NADH_copy, ignore_keys={field_name})
    assert equivalent
    assert comparison["same"] == set(NADH.__getstate__().keys()).difference(
        {field_name}
    )
    assert comparison["modified"] == {}
    assert comparison["added"] == set()
    assert comparison["removed"] == set()
    assert NADH.__getattribute__(field_name) != NADH_copy.__getattribute__(field_name)
    NADH.__setattr__(field_name, "Test")
    equivalent, comparison = compare_state(NADH, NADH_copy, ignore_keys={field_name})
    assert equivalent
    assert comparison["same"] == set(NADH.__getstate__().keys()).difference(
        {field_name}
    )
    assert comparison["modified"] == {}
    assert comparison["added"] == set()
    assert comparison["removed"] == set()
    assert NADH.__getattribute__(field_name) != NADH_copy.__getattribute__(field_name)
    NADH.__setattr__(field_name, "C21H27N7O14P")
    equivalent, comparison = compare_state(NADH, NADH_copy, ignore_keys={field_name})
    assert equivalent
    assert comparison["same"] == set(NADH.__getstate__().keys()).difference(
        {field_name}
    )
    assert comparison["modified"] == {}
    assert comparison["added"] == set()
    assert comparison["removed"] == set()
    assert NADH.__getattribute__(field_name) != NADH_copy.__getattribute__(field_name)
    NADH.__setattr__(field_name, "e")
    equivalent, comparison = compare_state(NADH, NADH_copy, ignore_keys={field_name})
    assert equivalent
    assert comparison["same"] == set(NADH.__getstate__().keys()).difference(
        {field_name}
    )
    assert comparison["modified"] == {}
    assert comparison["added"] == set()
    assert comparison["removed"] == set()
    assert NADH.__getattribute__(field_name) != NADH_copy.__getattribute__(field_name)


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
