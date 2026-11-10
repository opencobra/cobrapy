"""Test functionalities of FASTCC."""

from typing import Callable, List

import pytest

from cobra import Model, Reaction
from cobra.flux_analysis import fastcc, find_blocked_reactions


@pytest.fixture(scope="module")
def figure1_model() -> Model:
    """Generate a toy model as described in [1]_ figure 1.

    References
    ----------
    .. [1] Vlassis N, Pacheco MP, Sauter T (2014)
           Fast Reconstruction of Compact Context-Specific Metabolic Network
           Models.
           PLoS Comput Biol 10(1): e1003424. doi:10.1371/journal.pcbi.1003424

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


def test_fastcc_benchmark(
    model: Model, benchmark: Callable, all_solvers: List[str]
) -> None:
    """Benchmark fastcc."""
    model.solver = all_solvers
    benchmark(fastcc, model)


def test_figure1(figure1_model: Model, all_solvers: List[str]) -> None:
    """Test FASTCC."""
    figure1_model.solver = all_solvers
    consistent_model = fastcc(figure1_model)
    expected_reactions = {"v1", "v3", "v4", "v5", "v6"}
    assert expected_reactions == {rxn.id for rxn in consistent_model.reactions}


def test_opposing(opposing_model: Model, all_solvers: List[str]) -> None:
    """Test FASTCC."""
    opposing_model.solver = all_solvers
    consistent_model = fastcc(opposing_model)
    expected_reactions = {"v1", "v2", "v3", "v4"}
    assert expected_reactions == {rxn.id for rxn in consistent_model.reactions}


def test_fastcc_against_nonblocked_rxns(model: Model, all_solvers: List[str]) -> None:
    """Test non-blocked reactions obtained by FASTCC."""
    model.solver = all_solvers
    model.tolerance = 1e-6
    fastcc_consistent_model = fastcc(model, 1e-3, 1e-6)
    blocked = find_blocked_reactions(model)
    fastcc_ids = set(rxn.id for rxn in fastcc_consistent_model.reactions)
    assert len(model.reactions) - len(blocked) == len(fastcc_consistent_model.reactions)
    assert fastcc_ids & set(blocked) == set()


def test_fastcc_against_nonblocked_rxns_large(large_model: Model) -> None:
    """Test non-blocked reactions obtained by FASTCC."""
    model = large_model
    model.tolerance = 1e-7
    fastcc_consistent_model = fastcc(model, 1e-3, 1e-7)
    blocked = find_blocked_reactions(model, zero_cutoff=1e-7)
    fastcc_ids = set(rxn.id for rxn in fastcc_consistent_model.reactions)
    assert len(model.reactions) - len(blocked) == len(fastcc_consistent_model.reactions)
    assert fastcc_ids & set(blocked) == set()


@pytest.mark.filterwarnings("error:Solver status is")
@pytest.mark.parametrize("reverse_mask", range(8))
@pytest.mark.parametrize("minimum_flux,close_path", [(0, False), (0, True), (1, False)])
def test_fastcc_signed_path(
    reverse_mask: int, minimum_flux: int, close_path: bool, all_solvers: str
) -> None:
    """Changing reaction orientation preserves the exact consistency set."""
    test_model = Model("signed_path")
    uptake, path, drain, dead_end, closed = [
        Reaction(rid) for rid in ("uptake", "path", "drain", "dead_end", "closed")
    ]
    test_model.add_reactions([uptake, path, drain, dead_end, closed])
    uptake.reaction = "-> A"
    path.reaction = "A -> B"
    drain.reaction = "B ->"
    dead_end.reaction = "A -> C"
    closed.reaction = "A -> B"
    for reaction in test_model.reactions:
        reaction.bounds = (0, 10)
    closed.bounds = (0, 0)
    for reaction in (uptake, path, drain):
        reaction.bounds = (minimum_flux, 10)
    if close_path:
        path.bounds = (0, 0)
    for index, reaction in enumerate((uptake, path, drain)):
        if reverse_mask & (1 << index):
            reaction *= -1
    test_model.solver = all_solvers
    test_model.objective = drain
    before_objective = str(test_model.objective.expression)
    before = [(r.id, r.reaction, r.bounds) for r in test_model.reactions]
    before_variables = set(test_model.variables.keys())
    before_constraints = set(test_model.constraints.keys())

    consistent_model = fastcc(test_model)

    # C has no consumer, so dead_end is zero. The other mass balances force
    # equal physical pathway fluxes; all can reach 10 unless path is closed.
    expected = set() if close_path else {"uptake", "path", "drain"}
    assert {r.id for r in consistent_model.reactions} == expected
    assert [(r.id, r.reaction, r.bounds) for r in test_model.reactions] == before
    assert str(test_model.objective.expression) == before_objective
    assert set(test_model.variables.keys()) == before_variables
    assert set(test_model.constraints.keys()) == before_constraints


def test_fastcc_reversible_cycle_with_closed_reaction(all_solvers: str) -> None:
    """Reversible loops remain flux-consistent and fixed-zero reactions do not."""
    test_model = Model("reversible_cycle")
    forward, backward, closed = [Reaction(rid) for rid in ("f", "b", "closed")]
    test_model.add_reactions([forward, backward, closed])
    forward.reaction = "A <-> B"
    backward.reaction = "B <-> A"
    closed.reaction = "A -> B"
    forward.bounds = backward.bounds = (-10, 10)
    closed.bounds = (0, 0)
    test_model.solver = all_solvers
    assert {r.id for r in fastcc(test_model).reactions} == {"f", "b"}
