"""Test functionalities of MOMA."""

from typing import List

import numpy as np
import pytest

from cobra import Model
from cobra.flux_analysis.moma import add_moma


def test_moma_sanity(model: Model, qp_solvers: List[str]) -> None:
    """Test optimization criterion and optimality for MOMA."""
    model.solver = qp_solvers
    sol = model.optimize()

    with model:
        model.reactions.PFK.knock_out()
        knock_sol = model.optimize()
        ssq = (knock_sol.fluxes - sol.fluxes).pow(2).sum()

    with model:
        add_moma(model, linear=False)
        model.reactions.PFK.knock_out()
        moma_sol = model.optimize()
        moma_ssq = (moma_sol.fluxes - sol.fluxes).pow(2).sum()

    # Use normal FBA as reference solution.
    with model:
        add_moma(model, solution=sol, linear=False)
        model.reactions.PFK.knock_out()
        moma_ref_sol = model.optimize()
        moma_ref_ssq = (moma_ref_sol.fluxes - sol.fluxes).pow(2).sum()

    assert np.isclose(moma_sol.objective_value, moma_ssq)
    assert moma_ssq < ssq
    assert np.isclose(moma_sol.objective_value, moma_ref_sol.objective_value)
    assert np.isclose(moma_ssq, moma_ref_ssq)


def test_linear_moma_sanity(model: Model, all_solvers: List[str]) -> None:
    """Test optimization criterion and optimality for linear MOMA."""
    model.solver = all_solvers
    sol = model.optimize()

    with model:
        model.reactions.PFK.knock_out()
        knock_sol = model.optimize()
        sabs = (knock_sol.fluxes - sol.fluxes).abs().sum()

    with model:
        add_moma(model, linear=True)
        model.reactions.PFK.knock_out()
        moma_sol = model.optimize()
        moma_sabs = (moma_sol.fluxes - sol.fluxes).abs().sum()

    # Use normal FBA as reference solution.
    with model:
        add_moma(model, solution=sol, linear=True)
        model.reactions.PFK.knock_out()
        moma_ref_sol = model.optimize()
        moma_ref_sabs = (moma_ref_sol.fluxes - sol.fluxes).abs().sum()

    assert np.allclose(moma_sol.objective_value, moma_sabs)
    assert moma_sabs < sabs
    assert np.isclose(moma_sol.objective_value, moma_ref_sol.objective_value)
    assert np.isclose(moma_sabs, moma_ref_sabs)

    with model:
        add_moma(model, linear=True)
        with pytest.raises(ValueError):
            add_moma(model)
