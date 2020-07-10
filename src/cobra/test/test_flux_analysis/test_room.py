# -*- coding: utf-8 -*-

"""Test functionalities of ROOM."""

from __future__ import absolute_import

import numpy as np

from cobra.flux_analysis.room import add_room


def test_room_sanity(model, all_solvers):
    """Test optimization criterion and optimality for ROOM."""
    model.solver = all_solvers
    sol = model.optimize()
    with model:
        model.reactions.PYK.knock_out()
        knock_sol = model.optimize()

    with model:
        # Internally uses pFBA as reference solution.
        add_room(model)
        model.reactions.PYK.knock_out()
        room_sol = model.optimize()

    with model:
        # Use FBA as reference solution.
        add_room(model, solution=sol)
        model.reactions.PYK.knock_out()
        room_sol_ref = model.optimize()

    flux_change = (sol.fluxes - knock_sol.fluxes).abs().sum()
    flux_change_room = (sol.fluxes - room_sol.fluxes).abs().sum()
    flux_change_room_ref = (sol.fluxes - room_sol_ref.fluxes).abs().sum()
    # Expect the ROOM solution to have smaller flux changes in
    # reactions compared to a normal FBA.
    assert flux_change_room < flux_change or np.isclose(
        flux_change_room, flux_change, atol=1e-06
    )
    # Expect the FBA-based reference to have less change in
    # flux distribution.
    assert flux_change_room_ref > flux_change_room or np.isclose(
        flux_change_room_ref, flux_change_room, atol=1e-06
    )


def test_linear_room_sanity(model, all_solvers):
    """Test optimization criterion and optimality for linear ROOM."""
    model.solver = all_solvers
    sol = model.optimize()
    with model:
        model.reactions.PYK.knock_out()
        knock_sol = model.optimize()

    with model:
        # Internally uses pFBA as reference solution.
        add_room(model, linear=True)
        model.reactions.PYK.knock_out()
        room_sol = model.optimize()

    with model:
        # Use FBA as reference solution.
        add_room(model, solution=sol, linear=True)
        model.reactions.PYK.knock_out()
        room_sol_ref = model.optimize()

    flux_change = (sol.fluxes - knock_sol.fluxes).abs().sum()
    flux_change_room = (sol.fluxes - room_sol.fluxes).abs().sum()
    flux_change_room_ref = (sol.fluxes - room_sol_ref.fluxes).abs().sum()
    # Expect the ROOM solution to have smaller flux changes in
    # reactions compared to a normal FBA.
    assert flux_change_room < flux_change or np.isclose(
        flux_change_room, flux_change, atol=1e-06
    )
    # Expect the FBA-based reference to have less change in
    # flux distribution.
    assert flux_change_room_ref > flux_change_room or np.isclose(
        flux_change_room_ref, flux_change_room, atol=1e-06
    )
