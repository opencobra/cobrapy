"""Test functionalities of the SPECTRA consistency check."""

from typing import Callable, List

import pytest

from cobra import Model, Reaction
from cobra.flux_analysis import fastcc, find_blocked_reactions, spectra_cc
from cobra.flux_analysis.spectra_cc import _find_consistent_reaction_ids


@pytest.fixture(scope="module")
def figure1_model() -> Model:
    """Generate a toy model with a single blocked reaction.

    The reaction ``v2`` is the only producer and consumer of ``B``, so it
    cannot carry flux at steady state.

    """
    test_model = Model("figure 1")
    v1 = Reaction("v1")
    v2 = Reaction("v2")
    v3 = Reaction("v3")
    v4 = Reaction("v4")
    v5 = Reaction("v5")
    v6 = Reaction("v6")

    test_model.add_reactions([v1, v2, v3, v4, v5, v6])

    v1.reaction = "-> 2 A"
    v2.reaction = "A <-> B"
    v3.reaction = "A -> D"
    v4.reaction = "A -> C"
    v5.reaction = "C -> D"
    v6.reaction = "D ->"

    v1.bounds = (0.0, 3.0)
    v2.bounds = (-3.0, 3.0)
    v3.bounds = (0.0, 3.0)
    v4.bounds = (0.0, 3.0)
    v5.bounds = (0.0, 3.0)
    v6.bounds = (0.0, 3.0)

    test_model.objective = v6
    return test_model


@pytest.fixture(scope="module")
def opposing_model() -> Model:
    """Generate a toy model with opposing reversible reactions.

    This toy model ensures that two opposing reversible reactions do not
    appear as blocked.

    """
    test_model = Model("opposing")
    v1 = Reaction("v1")
    v2 = Reaction("v2")
    v3 = Reaction("v3")
    v4 = Reaction("v4")

    test_model.add_reactions([v1, v2, v3, v4])

    v1.reaction = "-> 2 A"
    v2.reaction = "A -> C"  # Later made reversible via bounds.
    v3.reaction = "D -> C"  # Later made reversible via bounds.
    v4.reaction = "D ->"

    v1.bounds = 0.0, 3.0
    v2.bounds = -3.0, 3.0
    v3.bounds = -3.0, 3.0
    v4.bounds = 0.0, 3.0

    test_model.objective = v4
    return test_model


@pytest.fixture(scope="module")
def backwards_model() -> Model:
    """Generate a toy model whose reactions only carry negative flux.

    Every reaction is restricted to a non-positive flux, which exercises the
    re-orientation SPECTRA applies in place of flipping the stoichiometric
    matrix. Read in reverse, ``v1`` to ``v3`` form a source-to-sink pathway
    over ``A`` and ``B``. ``v4`` is blocked because ``C`` has no consumer, and
    ``v5`` is blocked because it is closed outright.

    """
    test_model = Model("backwards")
    v1 = Reaction("v1")
    v2 = Reaction("v2")
    v3 = Reaction("v3")
    v4 = Reaction("v4")
    v5 = Reaction("v5")

    test_model.add_reactions([v1, v2, v3, v4, v5])

    v1.reaction = "A ->"  # Reads as "-> A" in reverse.
    v2.reaction = "B -> A"  # Reads as "A -> B" in reverse.
    v3.reaction = "-> B"  # Reads as "B ->" in reverse.
    v4.reaction = "C ->"  # Reads as "-> C" in reverse, and C has no consumer.
    v5.reaction = "D ->"

    v1.bounds = (-3.0, 0.0)
    v2.bounds = (-3.0, 0.0)
    v3.bounds = (-3.0, 0.0)
    v4.bounds = (-3.0, 0.0)
    v5.bounds = (0.0, 0.0)

    return test_model


def test_spectra_cc_benchmark(
    model: Model, benchmark: Callable, all_solvers: List[str]
) -> None:
    """Benchmark spectra_cc."""
    model.solver = all_solvers
    benchmark(spectra_cc, model, seed=0)


def test_figure1(figure1_model: Model, all_solvers: List[str]) -> None:
    """Test that the blocked reaction of the toy model is removed."""
    figure1_model.solver = all_solvers
    consistent_model = spectra_cc(figure1_model, seed=0)
    expected_reactions = {"v1", "v3", "v4", "v5", "v6"}
    assert expected_reactions == {rxn.id for rxn in consistent_model.reactions}


def test_opposing(opposing_model: Model, all_solvers: List[str]) -> None:
    """Test that opposing reversible reactions are not reported as blocked."""
    opposing_model.solver = all_solvers
    consistent_model = spectra_cc(opposing_model, seed=0)
    expected_reactions = {"v1", "v2", "v3", "v4"}
    assert expected_reactions == {rxn.id for rxn in consistent_model.reactions}


def test_backwards(backwards_model: Model, all_solvers: List[str]) -> None:
    """Test reactions restricted to negative flux."""
    backwards_model.solver = all_solvers
    consistent_model = spectra_cc(backwards_model, seed=0)
    expected_reactions = {"v1", "v2", "v3"}
    assert expected_reactions == {rxn.id for rxn in consistent_model.reactions}


def test_spectra_cc_against_nonblocked_rxns(
    model: Model, all_solvers: List[str]
) -> None:
    """Test non-blocked reactions obtained by SPECTRA."""
    model.solver = all_solvers
    model.tolerance = 1e-6
    spectra_consistent_model = spectra_cc(model, 1e-3, 1e-6, seed=0)
    blocked = find_blocked_reactions(model)
    spectra_ids = {rxn.id for rxn in spectra_consistent_model.reactions}
    assert len(model.reactions) - len(blocked) == len(
        spectra_consistent_model.reactions
    )
    assert spectra_ids & set(blocked) == set()


def test_spectra_cc_against_fastcc(model: Model, all_solvers: List[str]) -> None:
    """Test that SPECTRA and FASTCC agree on the consistent reactions."""
    model.solver = all_solvers
    model.tolerance = 1e-6
    spectra_consistent_model = spectra_cc(model, 1e-3, 1e-6, seed=0)
    fastcc_consistent_model = fastcc(model, 1e-3, 1e-6)
    assert {rxn.id for rxn in spectra_consistent_model.reactions} == {
        rxn.id for rxn in fastcc_consistent_model.reactions
    }


def test_spectra_cc_is_seeded(model: Model, all_solvers: List[str]) -> None:
    """Test that the randomized objective does not change the outcome."""
    model.solver = all_solvers
    model.tolerance = 1e-6
    first, _ = _find_consistent_reaction_ids(model, 1e-3, 1e-6, seed=0)
    repeated, _ = _find_consistent_reaction_ids(model, 1e-3, 1e-6, seed=0)
    other, _ = _find_consistent_reaction_ids(model, 1e-3, 1e-6, seed=7)
    assert first == repeated
    assert first == other


def test_spectra_cc_flux_threshold_semantics(model: Model) -> None:
    """Test that a cutoff at the flux threshold drops low-flux reactions.

    Passing ``zero_cutoff`` at the flux threshold reproduces the reference
    MATLAB implementation, which requires a reaction to reach the threshold
    itself rather than to merely carry flux. Reactions whose attainable flux
    falls between the two cutoffs are then dropped as well.

    """
    model.tolerance = 1e-7
    any_flux, _ = _find_consistent_reaction_ids(model, 1e-3, 1e-7, seed=0)
    at_threshold, _ = _find_consistent_reaction_ids(model, 1e-3, 9.9e-4, seed=0)
    assert at_threshold <= any_flux


def test_spectra_cc_rejects_cutoff_below_tolerance(model: Model) -> None:
    """Test that a cutoff below the solver tolerance is rejected."""
    model.tolerance = 1e-6
    with pytest.raises(ValueError):
        spectra_cc(model, 1e-3, 1e-9)


def test_spectra_cc_against_nonblocked_rxns_large(large_model: Model) -> None:
    """Test non-blocked reactions obtained by SPECTRA on a large model."""
    model = large_model
    model.tolerance = 1e-7
    spectra_consistent_model = spectra_cc(model, 1e-3, 1e-7, seed=0)
    blocked = find_blocked_reactions(model, zero_cutoff=1e-7)
    spectra_ids = {rxn.id for rxn in spectra_consistent_model.reactions}
    assert len(model.reactions) - len(blocked) == len(
        spectra_consistent_model.reactions
    )
    assert spectra_ids & set(blocked) == set()
