# -*- coding: utf-8 -*-

"""Test functionalities of MetaboliteSummary."""

from __future__ import absolute_import

import pytest

from cobra.flux_analysis.parsimonious import pfba
from cobra.test.test_core.test_summary import (
    captured_output, check_in_line, check_line)


@pytest.mark.parametrize("met", ["q8_c"])
def test_metabolite_summary_to_table_previous_solution(model, opt_solver, met):
    """Test metabolite summary._to_table() of previous solution."""
    model.solver = opt_solver
    solution = pfba(model)

    expected_entry = ['PRODUCING CYTBD     100    43.6  2.0 h_c + 0.5 o2_c + '
                      'q8h2_c --> h2o_c + 2.0 h_...']

    with captured_output() as (out, _):
        print(model.metabolites.get_by_id(met).summary(solution))
        check_in_line(out.getvalue(), expected_entry)


@pytest.mark.parametrize("met", ["q8_c"])
def test_metabolite_summary_to_frame_previous_solution(model, opt_solver, met):
    """Test metabolite summary.to_frame() of previous solution."""
    model.solver = opt_solver
    solution = pfba(model)

    expected_percent = [100.0, 88.4, 11.6]

    out_df = model.metabolites.get_by_id(met).summary(solution).to_frame()

    assert out_df['PERCENT'].round(1).tolist() == expected_percent


@pytest.mark.parametrize("met, names", [
    ("q8_c", False),
    ("q8_c", True)
])
def test_metabolite_summary_to_table(model, opt_solver, met, names):
    """Test metabolite summary._to_table()."""
    model.solver = opt_solver
    model.optimize()

    if names:
        expected_entry = ['PRODUCING cytochrome oxidase bd '
                          '(ubiquinol-8: 2 protons)    100    43.6  '
                          '2.0 H+ + 0.5 O2 + Ubiquinol-8 --> H2O + 2.0 H+...']
    else:
        expected_entry = ['PRODUCING CYTBD     100    43.6  2.0 h_c + '
                          '0.5 o2_c + q8h2_c --> h2o_c + 2.0 h_...']

    with captured_output() as (out, _):
        print(model.metabolites.get_by_id(met).summary(names=names))
        check_in_line(out.getvalue(), expected_entry)


@pytest.mark.parametrize("met, names", [
    ("q8_c", False),
    ("q8_c", True)
])
def test_metabolite_summary_to_frame(model, opt_solver, met, names):
    """Test metabolite summary.to_frame()."""
    model.solver = opt_solver
    model.optimize()

    expected_percent = [100.0, 88.4, 11.6]

    out_df = model.metabolites.get_by_id(met).summary(names=names).to_frame()

    assert out_df['PERCENT'].round(1).tolist() == expected_percent


@pytest.mark.parametrize("fraction, met", [(0.99, "fdp_c")])
def test_metabolite_summary_to_table_with_fva(model, opt_solver, fraction,
                                              met):
    """Test metabolite summary._to_table() (using FVA)."""
    #     pytest.xfail("FVA currently buggy")

    model.solver = opt_solver
    model.optimize()

    expected_entry = ['PRODUCING PFK   100     7.48    6.17      9.26    '
                      'atp_c + f6p_c --> adp_c + fdp_c + h_c']

    with captured_output() as (out, _):
        print(model.metabolites.get_by_id(met).summary(fva=fraction))
        check_line(out.getvalue(), expected_entry)


@pytest.mark.parametrize("fraction, met", [(0.99, "fdp_c")])
def test_metabolite_summary_to_frame_with_fva(model, opt_solver, fraction,
                                              met):
    """Test metabolite summary.to_frame() (using FVA)."""

    model.solver = opt_solver
    model.optimize()

    expected_percent = [100.0, 100.0, 0.0]

    out_df = model.metabolites.get_by_id(met).summary(fva=fraction).to_frame()

    assert out_df['PERCENT'].tolist() == expected_percent
