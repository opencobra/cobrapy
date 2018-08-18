# -*- coding: utf-8 -*-

"""Test functionalities of pFBA."""

import warnings

import numpy as np

import pytest
from cobra.exceptions import Infeasible
from cobra.flux_analysis.parsimonious import add_pfba, pfba


def test_pfba_benchmark(large_model, benchmark, all_solvers):
    """Benchmark pFBA functionality."""
    large_model.solver = all_solvers
    benchmark(pfba, large_model)


def test_pfba(model, all_solvers):
    """Test pFBA functionality."""
    model.solver = all_solvers
    with model:
        add_pfba(model)
        with pytest.raises(ValueError):
            add_pfba(model)

    expression = model.objective.expression
    n_constraints = len(model.constraints)
    solution = pfba(model)
    assert solution.status == "optimal"
    assert np.isclose(solution.x_dict["Biomass_Ecoli_core"],
                      0.8739, atol=1e-4, rtol=0.0)
    abs_x = [abs(i) for i in solution.x]
    assert np.isclose(sum(abs_x), 518.4221, atol=1e-4, rtol=0.0)
    # test changes to model reverted
    assert expression == model.objective.expression
    assert len(model.constraints) == n_constraints

    # needed?
    # Test desired_objective_value
    # desired_objective = 0.8
    # pfba(model, solver=solver,
    #                       desired_objective_value=desired_objective)
    # abs_x = [abs(i) for i in model.solution.x]
    # assert model.solution.status == "optimal"
    # assert abs(model.solution.f - desired_objective) < 0.001
    # assert abs(sum(abs_x) - 476.1594) < 0.001

    # TODO: parametrize fraction (DRY it up)
    # Test fraction_of_optimum
    solution = pfba(model, fraction_of_optimum=0.95)
    assert solution.status == "optimal"
    assert np.isclose(solution.x_dict["Biomass_Ecoli_core"],
                      0.95 * 0.8739, atol=1e-4, rtol=0.0)
    abs_x = [abs(i) for i in solution.x]
    assert np.isclose(sum(abs_x), 493.4400, atol=1e-4, rtol=0.0)

    # Infeasible solution
    model.reactions.ATPM.lower_bound = 500
    with warnings.catch_warnings():
        warnings.simplefilter("error", UserWarning)
        with pytest.raises((UserWarning, Infeasible, ValueError)):
            pfba(model)
