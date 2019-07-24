# -*- coding: utf-8 -*-

"""Test functionalities of ModelSummary."""

from __future__ import absolute_import

import numpy as np
import pytest

from cobra.test.test_core.test_summary import (
    captured_output, check_in_line, check_line)


@pytest.mark.parametrize("names", [False, True])
def test_model_summary_to_table_previous_solution(model, opt_solver, names):
    """Test Summary._to_table() of previous solution."""
    model.solver = opt_solver
    solution = model.optimize()
    rxn_test = model.exchanges[0]

    if names:
        met_test = list(rxn_test.metabolites.keys())[0].name
    else:
        met_test = list(rxn_test.metabolites.keys())[0].id

    solution.fluxes[rxn_test.id] = 321

    with captured_output() as (out, _):
        print(model.summary(solution, names=names))
        check_in_line(out.getvalue(), [met_test + '321'])


@pytest.mark.parametrize("names", [False, True])
def test_model_summary_to_frame_previous_solution(model, opt_solver, names):
    """Test Summary.to_frame() of previous solution."""
    model.solver = opt_solver
    solution = model.optimize()
    rxn_test = model.exchanges[0]
    solution.fluxes[rxn_test.id] = 321
    out_df = model.summary(solution, names=names).to_frame()

    assert out_df.loc[0, ('OUT_FLUXES', 'FLUX')] == 321


@pytest.mark.parametrize("names", [False, True])
def test_model_summary_to_table(model, opt_solver, names):
    """Test model.summary()._to_table()."""
    model.solver = opt_solver
    # test non-fva version (these should be fixed for textbook model)
    if names:
        expected_entry = ['        O2    21.8       H2O       29.2     '
                          'Biomass Objective Function with GAM   0.874   ']
    else:
        expected_entry = ['     o2_e    21.8      h2o_e      29.2     '
                          'Biomass_Ecoli_core   0.874   ']

    model.optimize()

    with captured_output() as (out, _):
        print(model.summary(names=names))
        check_in_line(out.getvalue(), expected_entry)


@pytest.mark.parametrize("names", [False, True])
def test_model_summary_to_frame(model, opt_solver, names):
    """Test model.summary().to_frame()."""
    model.solver = opt_solver
    # test non-fva version (these should be fixed for textbook model)
    if names:
        expected_in_fluxes = ['O2', 'D-Glucose', 'Ammonium', 'Phosphate']
        expected_out_fluxes = ['H2O', 'CO2', 'H+', np.nan]
    else:
        expected_in_fluxes = ['o2_e', 'glc__D_e', 'nh4_e', 'pi_e']
        expected_out_fluxes = ['h2o_e', 'co2_e', 'h_e', np.nan]

    model.optimize()

    out_df = model.summary(names=names).to_frame()

    assert out_df[('IN_FLUXES', 'ID')].tolist() == expected_in_fluxes
    assert out_df[('OUT_FLUXES', 'ID')].tolist() == expected_out_fluxes


@pytest.mark.parametrize("fraction", [0.95])
def test_model_summary_to_table_with_fva(model, opt_solver, fraction):
    """Test model summary._to_table() (using FVA)."""
    if opt_solver == "optlang-gurobi":
        pytest.xfail("FVA currently buggy")
    # test non-fva version (these should be fixed for textbook model)
    expected_entry = ['     o2_e    21.8      19.9      23.7       h2o_e     '
                      '29.2         25       30.7     Biomass_Ecoli_core   '
                      '0.874   ']

    model.solver = opt_solver
    solution = model.optimize()

    with captured_output() as (out, _):
        print(model.summary(solution, fva=fraction))
        check_in_line(out.getvalue(), expected_entry)


@pytest.mark.parametrize("fraction", [0.95])
def test_model_summary_to_frame_with_fva(model, opt_solver, fraction):
    """Test model summary.to_frame() (using FVA)."""
    if opt_solver == "optlang-gurobi":
        pytest.xfail("FVA currently buggy")
    # test non-fva version (these should be fixed for textbook model)
    expected_in_fluxes = ['o2_e', 'glc__D_e', 'nh4_e', 'pi_e', np.nan,
                          np.nan, np.nan, np.nan, np.nan, np.nan, np.nan,
                          np.nan]
    expected_out_fluxes = ['h2o_e', 'co2_e', 'h_e', 'for_e', 'ac_e', 'acald_e',
                           'pyr_e', 'etoh_e', 'lac__D_e', 'succ_e', 'akg_e',
                           'glu__L_e']

    model.solver = opt_solver
    solution = model.optimize()
    out_df = model.summary(solution, fva=fraction).to_frame()

    assert out_df[('IN_FLUXES', 'ID')].tolist() == expected_in_fluxes
    assert out_df[('OUT_FLUXES', 'ID')].tolist() == expected_out_fluxes
