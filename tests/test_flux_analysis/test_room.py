"""Test functionalities of ROOM."""

from typing import List

import pytest

from cobra import Model
from cobra.flux_analysis.room import add_room


@pytest.mark.parametrize(
    "linear, delta, eps",
    [(True, 0.03, 1e-3), (False, 0.03, 1e-3), (True, 0.1, 1e-2), (False, 0.1, 1e-2)],
)
def test_room_sanity(
    model: Model, all_solvers: List[str], linear: bool, delta: float, eps: float
) -> None:
    """Test optimization criterion and optimality for ROOM."""
    model.solver = all_solvers
    sol = model.optimize()
    with model:
        model.reactions.PYK.knock_out()
        knock_sol = model.optimize()

    with model:
        # let it calculate its own reference solution (pFBA)
        add_room(model, linear=linear, delta=delta, epsilon=eps)
        model.reactions.PYK.knock_out()
        room_sol = model.optimize()

    with model:
        # use the more distant FBA reference
        add_room(model, solution=sol, linear=linear, delta=delta, epsilon=eps)
        model.reactions.PYK.knock_out()
        room_sol_ref = model.optimize()

    # The
    flux_change = (sol.fluxes - knock_sol.fluxes).abs()
    flux_change_room = (sol.fluxes - room_sol.fluxes).abs()
    flux_change_room_ref = (sol.fluxes - room_sol_ref.fluxes).abs()
    rxn_count_naive = (flux_change > delta * sol.fluxes.abs() + eps + 1e-6).sum()
    rxn_count_room = (flux_change_room > delta * sol.fluxes.abs() + eps + 1e-6).sum()
    rxn_count_room_ref = (
        flux_change_room_ref > delta * sol.fluxes.abs() + eps + 1e-6
    ).sum()

    # Expect the ROOM solution to have less changed reactions then a pFBA or FBA
    assert rxn_count_room <= rxn_count_naive
    assert rxn_count_room_ref <= rxn_count_naive
