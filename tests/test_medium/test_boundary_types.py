"""Test functionalities of boundary type detection functions."""


import pytest

from cobra.core import Metabolite, Model, Reaction
from cobra.medium import (
    find_boundary_types,
    find_external_compartment,
    is_boundary_type,
)


def test_find_external_compartment_single(model: Model) -> None:
    """Test external compartment identification."""
    # by name
    assert find_external_compartment(model) == "e"
    # from boundary counts
    for m in model.metabolites:
        if m.compartment == "e":
            m.compartment = "outside"
    for r in model.reactions:
        r._compartments = None
        assert find_external_compartment(model) == "outside"
    # names are always right
    model.exchanges[0].reactants[0].compartment = "extracellular"
    assert find_external_compartment(model) == "extracellular"


def test_find_external_compartment_multi(model: Model) -> None:
    """Test multiple external compartment identification."""
    for r in model.reactions:
        r._compartments = None
    model.exchanges[0].reactants[0].compartment = "extracellular"
    # still works due to different boundary numbers
    assert find_external_compartment(model) == "e"
    model.exchanges[1].reactants[0].compartment = "extra cellular"
    model.remove_reactions(model.exchanges)
    # Now fails because same boundary count
    with pytest.raises(RuntimeError):
        find_external_compartment(model)


def test_no_names_or_boundary_reactions(empty_model: Model) -> None:
    """Test absence of name or boundary reactions."""
    with pytest.raises(RuntimeError):
        find_external_compartment(empty_model)


def test_find_boundary_types_exchange(model: Model) -> None:
    """Test boundary type identification for exchanges."""
    ex = model.exchanges
    assert all(r.id.startswith("EX_") for r in ex)
    ex = find_boundary_types(model, "exchange", "e")
    assert all(r.id.startswith("EX_") for r in ex)


def test_find_boundary_types_demand(model: Model) -> None:
    """Test boundary type identification for demands."""
    dm = Reaction("demand")
    model.add_reactions([dm])
    dm.build_reaction_from_string("atp_c ->")
    dm = model.demands
    assert len(dm) == 1
    assert "demand" in [r.id for r in dm]


def test_find_boundary_types_sink(model: Model) -> None:
    """Test boundary type identification for sinks."""
    sn = Reaction("sink")
    model.add_reactions([sn])
    sn.build_reaction_from_string("atp_c <->")
    sn.bounds = -1000, 1000
    sn = model.sinks
    assert len(sn) == 1
    assert "sink" in [r.id for r in sn]


def test_no_boundary_reactions(empty_model: Model) -> None:
    """Test proper identification of no boundary reactions."""
    assert find_boundary_types(empty_model, "e", None) == []


def test_is_boundary_type(model: Model) -> None:
    """Test correct identification of boundary types for reactions."""
    assert not is_boundary_type(model.reactions.ATPM, "exchange", "e")
    model.reactions.ATPM.annotation["sbo"] = "SBO:0000627"
    assert is_boundary_type(model.reactions.ATPM, "exchange", "bla")
    model.reactions.ATPM.annotation["sbo"] = "SBO:0000632"
    assert not is_boundary_type(model.reactions.ATPM, "exchange", "e")


def test_bad_exchange(model: Model) -> None:
    """Test bad exchange reaction identification."""
    with pytest.raises(ValueError):
        m = Metabolite("baddy", compartment="nonsense")
        model.add_boundary(m, type="exchange")
    m = Metabolite("goody", compartment="e")
    rxn = model.add_boundary(m, type="exchange")
    assert isinstance(rxn, Reaction)
