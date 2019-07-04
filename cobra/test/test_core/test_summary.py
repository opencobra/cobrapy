# -*- coding: utf-8 -*-

"""Test functionalities of summary methods."""

from __future__ import absolute_import

import re
import sys
from contextlib import contextmanager

import numpy as np
import pytest
from six import StringIO

from cobra.flux_analysis.parsimonious import pfba


# Helper functions for the module

@contextmanager
def captured_output():
    """A context manager to test the IO summary methods."""
    new_out, new_err = StringIO(), StringIO()
    old_out, old_err = sys.stdout, sys.stderr
    try:
        sys.stdout, sys.stderr = new_out, new_err
        yield sys.stdout, sys.stderr

    finally:
        sys.stdout, sys.stderr = old_out, old_err


def check_line(output, expected_entries,
               pattern=re.compile(r"\s")):
    """Ensure each expected entry is in the output."""
    output_set = set(
        pattern.sub("", line) for line in output.splitlines())
    for elem in expected_entries:
        assert pattern.sub("", elem) in output_set


def check_in_line(output, expected_entries,
                  pattern=re.compile(r"\s")):
    """Ensure each expected entry is contained in the output."""
    output_strip = [pattern.sub("", line) for line in
                    output.splitlines()]
    for elem in expected_entries:
        assert any(
            pattern.sub("", elem) in line for line in output_strip), \
            "Not found: {} in:\n{}".format(pattern.sub("", elem),
                                           "\n".join(output_strip))


# Test functions

@pytest.mark.parametrize("names", [False, True])
def test_model_summary_to_table_previous_solution(model, opt_solver, names):
    """Test Summary.to_table() of previous solution."""
    model.solver = opt_solver
    solution = model.optimize()
    rxn_test = model.exchanges[0]
    if names:
        met_test = list(rxn_test.metabolites.keys())[0].name
    else:
        met_test = list(rxn_test.metabolites.keys())[0].id

    solution.fluxes[rxn_test.id] = 321

    with captured_output() as (out, _):
        model.summary(solution, names=names).to_table()
    check_in_line(out.getvalue(), [met_test + '321'])


@pytest.mark.parametrize("names", [False, True])
def test_model_summary_to_frame_previous_solution(model, opt_solver, names):
    """Test Summary.to_frame() of previous solution."""
    model.solver = opt_solver
    solution = model.optimize()
    rxn_test = model.exchanges[0]
    if names:
        met_test = list(rxn_test.metabolites.keys())[0].name
    else:
        met_test = list(rxn_test.metabolites.keys())[0].id

    solution.fluxes[rxn_test.id] = 321

    out_df = model.summary(solution, names=names).to_frame()

    assert out_df.loc[0, ('OUT_FLUXES', 'FLUX')] == 321


@pytest.mark.parametrize("names", [False, True])
def test_model_summary_to_table(model, opt_solver, names):
    """Test model.summary().to_table()."""
    model.solver = opt_solver
    # test non-fva version (these should be fixed for textbook model)
    if names:
        expected_entries = [
            'O2      21.8',
            'D-Glucose  10',
            'Ammonium      4.77',
            'Phosphate       3.21',
            'H2O  29.2',
            'CO2  22.8',
            'H+    17.5',
            'Biomass_Ecol...  0.874',
        ]
    else:
        expected_entries = [
            'o2_e      21.8',
            'glc__D_e  10',
            'nh4_e      4.77',
            'pi_e       3.21',
            'h2o_e  29.2',
            'co2_e  22.8',
            'h_e    17.5',
            'Biomass_Ecol...  0.874',
        ]
    # Need to use a different method here because
    # there are multiple entries per line.
    model.optimize()
    with captured_output() as (out, _):
        model.summary(names=names).to_table()
    check_in_line(out.getvalue(), expected_entries)

    # with model:
    #     model.objective = model.exchanges[0]
    #     model.summary()


@pytest.mark.parametrize("names", [False, True])
def test_model_summary_to_frame(model, opt_solver, names):
    """Test model.summary().to_frame()."""
    model.solver = opt_solver
    # test non-fva version (these should be fixed for textbook model)
    if names:
        expected_in_fluxes = ['O2', 'D-Glucose', 'Ammonium', 'Phosphate']
        expected_out_fluxes = ['H2O', 'CO2', 'H+', '']
    else:
        expected_in_fluxes = ['o2_e', 'glc__D_e', 'nh4_e', 'pi_e']
        expected_out_fluxes = ['h2o_e', 'co2_e', 'h_e', '']

    model.optimize()

    out_df = model.summary(names=names).to_frame()

    assert out_df[('IN_FLUXES', 'ID')].tolist() == expected_in_fluxes
    assert out_df[('OUT_FLUXES', 'ID')].tolist() == expected_out_fluxes


@pytest.mark.parametrize("fraction", [0.95])
def test_model_summary_to_table_with_fva(model, opt_solver, fraction):
    """Test model summary.to_table() (using FVA)."""
    if opt_solver == "optlang-gurobi":
        pytest.xfail("FVA currently buggy")
    # test non-fva version (these should be fixed for textbook model)
    expected_entries = [
        'idFluxRangeidFluxRangeBiomass_Ecol...0.874',
        'o2_e       21.8   [19.9, 23.7]'
        'h2o_e       29.2  [25, 30.7]',
        'glc__D_e   10     [9.52, 10]'
        'co2_e       22.8  [18.9, 24.7]',
        'nh4_e       4.77  [4.53, 5.16]'
        'h_e         17.5  [16.7, 22.4]',
        'pi_e        3.21  [3.05, 3.21]'
        'for_e        0    [0, 5.72]',
        'ac_e         0    [0, 1.91]',
        'pyr_e        0    [0, 1.27]',
        'lac__D_e     0    [0, 1.07]',
        'succ_e       0    [0, 0.837]',
        'glu__L_e     0    [0, 0.636]',
        'akg_e        0    [0, 0.715]',
        'etoh_e       0    [0, 1.11]',
        'acald_e      0    [0, 1.27]',
    ]
    # Need to use a different method here because
    # there are multiple entries per line.
    model.solver = opt_solver
    solution = model.optimize()
    with captured_output() as (out, _):
        model.summary(solution, fva=fraction).to_table()
    check_in_line(out.getvalue(), expected_entries)


@pytest.mark.parametrize("fraction", [0.95])
def test_model_summary_to_frame_with_fva(model, opt_solver, fraction):
    """Test model summary.to_frame() (using FVA)."""
    if opt_solver == "optlang-gurobi":
        pytest.xfail("FVA currently buggy")
    # test non-fva version (these should be fixed for textbook model)
    expected_in_fluxes = ['o2_e', 'glc__D_e', 'nh4_e', 'pi_e', '', '', '', '',
                          '', '', '', '']
    expected_out_fluxes = ['h2o_e', 'co2_e', 'h_e', 'for_e', 'ac_e', 'acald_e',
                           'pyr_e', 'etoh_e', 'lac__D_e', 'succ_e', 'akg_e',
                           'glu__L_e']

    model.solver = opt_solver
    solution = model.optimize()
    out_df = model.summary(solution, fva=fraction).to_frame()

    assert out_df[('IN_FLUXES', 'ID')].tolist() == expected_in_fluxes
    assert out_df[('OUT_FLUXES', 'ID')].tolist() == expected_out_fluxes


@pytest.mark.parametrize("met", ["q8_c"])
def test_metabolite_summary_to_table_previous_solution(
        model, opt_solver, met):
    """Test metabolite summary.to_table() of previous solution."""
    model.solver = opt_solver
    solution = pfba(model)

    expected_entries = [
        'PRODUCING REACTIONS -- Ubiquinone-8 (q8_c)',
        '%       FLUX  RXN ID    REACTION',
        '100%   43.6   CYTBD     2.0 h_c + 0.5 o2_c + q8h2_c --> h2o_c + 2.0 '
        'h_e...',
        'CONSUMING REACTIONS -- Ubiquinone-8 (q8_c)',
        '%       FLUX  RXN ID    REACTION',
        '88%    38.5   NADH16    4.0 h_c + nadh_c + q8_c --> 3.0 h_e + nad_c +'
        ' q...',
        '12%     5.06  SUCDi     q8_c + succ_c --> fum_c + q8h2_c'
    ]

    with captured_output() as (out, _):
        model.metabolites.get_by_id(met).summary(solution).to_table()
    check_in_line(out.getvalue(), expected_entries)


@pytest.mark.parametrize("met", ["q8_c"])
def test_metabolite_summary_to_frame_previous_solution(
        model, opt_solver, met):
    """Test metabolite summary.to_frame() of previous solution."""
    model.solver = opt_solver
    solution = pfba(model)

    expected_percent = ['100%', '88%', '12%']

    out_df = model.metabolites.get_by_id(met).summary(solution).to_frame()

    assert out_df['PERCENT'].tolist() == expected_percent


@pytest.mark.parametrize("met, names", [
    ("q8_c", False),
    ("q8_c", True)
])
def test_metabolite_summary_to_table(model, opt_solver, met, names):
    """Test metabolite summary.to_table()."""
    model.solver = opt_solver
    model.optimize()
    with captured_output() as (out, _):
        model.metabolites.get_by_id(met).summary(names=names).to_table()

    if names:
        expected_entries = [
            'PRODUCING REACTIONS -- Ubiquinone-8 (q8_c)',
            '%       FLUX  RXN ID      REACTION',
            '100%   43.6   cytochr...  2.0 H+ + 0.5 O2 + Ubiquinol-8 --> '
            'H2O + 2.0 H+ ...',
            'CONSUMING REACTIONS -- Ubiquinone-8 (q8_c)',
            '%       FLUX  RXN ID      REACTION',
            '88%    38.5   NADH de...  4.0 H+ + Nicotinamide adenine '
            'dinucleotide - re...',
            '12%     5.06  succina...  Ubiquinone-8 + Succinate --> '
            'Fumarate + Ubiquin...'
        ]
    else:
        expected_entries = [
            'PRODUCING REACTIONS -- Ubiquinone-8 (q8_c)',
            '%       FLUX  RXN ID    REACTION',
            '100%   43.6   CYTBD     '
            '2.0 h_c + 0.5 o2_c + q8h2_c --> h2o_c + 2.0 h_e...',
            'CONSUMING REACTIONS -- Ubiquinone-8 (q8_c)',
            '%       FLUX  RXN ID    REACTION',
            '88%    38.5   NADH16    '
            '4.0 h_c + nadh_c + q8_c --> 3.0 h_e + nad_c + q...',
            '12%     5.06  SUCDi     q8_c + succ_c --> fum_c + q8h2_c',
        ]

    check_in_line(out.getvalue(), expected_entries)


@pytest.mark.parametrize("met, names", [
    ("q8_c", False),
    ("q8_c", True)
])
def test_metabolite_summary_to_frame(model, opt_solver, met, names):
    """Test metabolite summary.to_frame()."""
    model.solver = opt_solver
    model.optimize()

    expected_percent = ['100%', '88%', '12%']

    out_df = model.metabolites.get_by_id(met).summary(names=names).to_frame()

    assert out_df['PERCENT'].tolist() == expected_percent


@pytest.mark.parametrize("fraction, met", [(0.99, "fdp_c")])
def test_metabolite_summary_to_table_with_fva(model, opt_solver, fraction,
                                              met):
    """Test metabolite summary.to_table() (using FVA)."""
    #     pytest.xfail("FVA currently buggy")

    model.solver = opt_solver
    model.optimize()
    with captured_output() as (out, _):
        model.metabolites.get_by_id(met).summary(fva=fraction).to_table()

    expected_entries = [
        'PRODUCING REACTIONS -- D-Fructose 1,6-bisphosphate (fdp_c)',
        '%       FLUX  RANGE         RXN ID    REACTION',
        '100%    7.48  [6.17, 9.26]  PFK       '
        'atp_c + f6p_c --> adp_c + fdp_c + h_c',
        'CONSUMING REACTIONS -- D-Fructose 1,6-bisphosphate (fdp_c)',
        '%       FLUX  RANGE         RXN ID    REACTION',
        '100%    7.48  [6.17, 8.92]  FBA       fdp_c <=> dhap_c + g3p_c',
        '0%      0     [0, 1.72]     FBP       '
        'fdp_c + h2o_c --> f6p_c + pi_c',
    ]

    check_line(out.getvalue(), expected_entries)


@pytest.mark.parametrize("fraction, met", [(0.99, "fdp_c")])
def test_metabolite_summary_to_frame_with_fva(model, opt_solver, fraction,
                                              met):
    """Test metabolite summary.to_frame() (using FVA)."""

    model.solver = opt_solver
    model.optimize()

    expected_percent = ['100%', '100%', '0%']

    out_df = model.metabolites.get_by_id(met).summary(fva=fraction).to_frame()

    assert out_df['PERCENT'].tolist() == expected_percent


@pytest.mark.parametrize("rxn, names", [("ACALD", False), ("ACALD", True)])
def test_reaction_summary_to_table(model, rxn, names):
    """Test reaction summary.to_table()."""

    with captured_output() as (out, _):
        model.reactions.get_by_id(rxn).summary(names=names).to_table()

    if names:
        expected_entries = [
            'REACTION: Acetaldehyde + Coenzyme A + Nicotinamide adenine '
            'dinucleotide <=> Acetyl-CoA + H+ + Nicotinamide adenine '
            'dinucleotide - reduced',
            'GENES -- acetaldehyde dehydrogenase (acetylating) ACALD',
            'ID',
            'adhE',
            'mhpF',
            'REACTANTS -- acetaldehyde dehydrogenase (acetylating) ACALD',
            'ID                                   STOICHIOMETRIC COEFFICIENT'
            'COMPARTMENT',
            'Acetaldehyde                                                 -1  '
            'c',
            'Coenzyme A                                                   -1  '
            'c',
            'Nicotinamide adenine dinucleotide                            -1  '
            'c',
            'PRODUCTS -- acetaldehyde dehydrogenase (acetylating) ACALD',
            'ID                                             STOICHIOMETRIC '
            'COEFFICIENT  COMPARTMENT',
            'Acetyl-CoA                                                       '
            '       1  c',
            'H+                                                               '
            '       1  c',
            'Nicotinamide adenine dinucleotide - reduced                      '
            '       1  c'
        ]
    else:
        expected_entries = [
            'REACTION: acald_c + coa_c + nad_c <=> accoa_c + h_c + nadh_c',
            'GENES -- acetaldehyde dehydrogenase (acetylating) ACALD',
            'ID',
            'b1241',
            'b0351',
            'REACTANTS -- acetaldehyde dehydrogenase (acetylating) ACALD',
            'ID         STOICHIOMETRIC COEFFICIENT  COMPARTMENT',
            'acald_c                            -1  c',
            'coa_c                              -1  c',
            'nad_c                              -1  c',
            'PRODUCTS -- acetaldehyde dehydrogenase (acetylating) ACALD',
            'ID         STOICHIOMETRIC COEFFICIENT  COMPARTMENT',
            'accoa_c                             1  c',
            'h_c                                 1  c',
            'nadh_c                              1  c'
        ]

    check_line(out.getvalue(), expected_entries)


@pytest.mark.parametrize("rxn, names", [("ACALD", False), ("FUM", True)])
def test_reaction_summary_to_frame(model, rxn, names):
    """Test reaction summary.to_frame()."""

    out_df = model.reactions.get_by_id(rxn).summary(names=names).to_frame()

    if names:
        expected_gene_names = ['fumC', 'fumA', 'fumB']
        expected_met_names = ['Fumarate', 'H2O', 'L-Malate']

    else:
        expected_gene_names = ['b1241', 'b0351']
        expected_met_names = ['acald_c', 'coa_c', 'nad_c', 'accoa_c', 'h_c',
                              'nadh_c']

    assert all(out_df['REACTION', 'GENES', 'ID']
               .replace('', np.nan).dropna().tolist()) == \
        all(expected_gene_names)

    assert all(out_df['REACTION', 'METABOLITES', 'ID'].tolist()) == \
        all(expected_met_names)
