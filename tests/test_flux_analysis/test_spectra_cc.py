"""Test functionalities of the SPECTRA consistency check."""

from typing import Callable, List

import numpy as np
import pytest

from cobra import Model, Reaction
from cobra.flux_analysis import fastcc, find_blocked_reactions, spectra_cc
from cobra.flux_analysis.spectra_cc import (
    _find_consistent_reaction_ids,
    _find_flux_mode,
    _reaction_sign,
)


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


def test_spectra_cc_leaves_the_input_model_alone(model: Model) -> None:
    """Test that the auxiliary variables do not outlive the check."""
    model.tolerance = 1e-9
    before = (
        len(model.variables),
        len(model.constraints),
        len(model.reactions),
        str(model.objective.expression),
        model.objective.direction,
    )
    spectra_cc(model, 1e-3, 9.9e-4, seed=0)
    assert before == (
        len(model.variables),
        len(model.constraints),
        len(model.reactions),
        str(model.objective.expression),
        model.objective.direction,
    )


def test_spectra_cc_weights_ignore_iteration_order(
    model: Model, monkeypatch: pytest.MonkeyPatch
) -> None:
    """Test that each reaction gets the same weight however it was ordered.

    The objective coefficients are drawn in one block and zipped against the
    reactions positionally, so the order of the identifiers decides which
    reaction gets which weight. Callers pass a set, whose iteration order
    varies with ``PYTHONHASHSEED``, so a seeded run would otherwise differ
    between processes. The consistent set itself cannot show this, because
    the method is built to reach the same answer whatever the weights are.

    """
    model.tolerance = 1e-9
    rxn_ids = [rxn.id for rxn in model.reactions]
    signs = {rxn.id: _reaction_sign(rxn) for rxn in model.reactions}
    captured = []

    objective_class = model.problem.Objective
    original = objective_class.set_linear_coefficients

    def spy(self, coefficients, **kwargs):
        captured.append({var.name: value for var, value in coefficients.items()})
        return original(self, coefficients, **kwargs)

    monkeypatch.setattr(objective_class, "set_linear_coefficients", spy)

    for order in (rxn_ids, list(reversed(rxn_ids))):
        _find_flux_mode(model, order, signs, 1e-3, 9.9e-4, np.random.default_rng(0))

    assert len(captured) == 2
    assert len(captured[0]) == len(rxn_ids)
    assert captured[0] == captured[1]


def test_spectra_cc_keeps_every_attribute(large_model: Model) -> None:
    """Test that the returned model is the source, minus the blocked set.

    Building the result from the kept reactions rather than copying the whole
    model and deleting from it is what makes the check fast, so everything
    the source carried has to be shown to arrive: annotations on reactions,
    metabolites and genes, groups, the objective, and the element order.

    """
    large_model.tolerance = 1e-9
    consistent = spectra_cc(large_model, 1e-3, 9.9e-4, seed=0)
    kept = {rxn.id for rxn in consistent.reactions}
    assert 0 < len(kept) < len(large_model.reactions)

    for rxn in consistent.reactions:
        source = large_model.reactions.get_by_id(rxn.id)
        assert rxn.annotation == source.annotation
        assert rxn.notes == source.notes
        assert rxn.name == source.name
        assert rxn.subsystem == source.subsystem
        assert rxn.bounds == source.bounds
        assert rxn.gene_reaction_rule == source.gene_reaction_rule

    for met in consistent.metabolites:
        source = large_model.metabolites.get_by_id(met.id)
        assert met.annotation == source.annotation
        assert (met.name, met.formula, met.charge) == (
            source.name,
            source.formula,
            source.charge,
        )

    # The GPR is serialized as a string, so the genes are rebuilt by parsing
    # it and arrive carrying nothing but an identifier unless put back.
    annotated = [gene for gene in consistent.genes if gene.annotation]
    assert len(annotated) > 0
    for gene in consistent.genes:
        source = large_model.genes.get_by_id(gene.id)
        assert gene.annotation == source.annotation
        assert gene.name == source.name

    assert {group.id for group in consistent.groups} == {
        group.id for group in large_model.groups
    }
    assert len(large_model.groups) > 0
    for group in consistent.groups:
        source = large_model.groups.get_by_id(group.id)
        assert (group.name, group.kind) == (source.name, source.kind)
        assert {member.id for member in group.members} <= {
            member.id for member in source.members
        }

    assert consistent.objective.direction == large_model.objective.direction
    assert {
        rxn.id: rxn.objective_coefficient
        for rxn in consistent.reactions
        if rxn.objective_coefficient
    } == {
        rxn.id: rxn.objective_coefficient
        for rxn in large_model.reactions
        if rxn.objective_coefficient and rxn.id in kept
    }

    # Order follows the source, not the order the reactions introduced them.
    for collection in ("reactions", "metabolites", "genes"):
        produced = [obj.id for obj in getattr(consistent, collection)]
        expected = [
            obj.id
            for obj in getattr(large_model, collection)
            if obj.id in set(produced)
        ]
        assert produced == expected
