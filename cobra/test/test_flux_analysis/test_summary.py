# -*- coding: utf-8 -*-

"""Test functionalities of summary methods."""

from __future__ import absolute_import

import re
import sys
from contextlib import contextmanager

import pytest
from six import StringIO

from cobra.flux_analysis.parsimonious import pfba


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
        met_test = list(rxn_test.metabolites.keys())[0].nameop
    else:
        met_test = list(rxn_test.metabolites.keys())[0].id

    solution.fluxes[rxn_test.id] = 321

    out_df = model.summary(solution, names=names).to_frame()

    column_names = [['IN_FLUXES', 'OUT_FLUXES', 'OBJECTIVES'], ['ID', 'FLUX']]

    if names:
        data = np.array([['O2', 21.799492655998794, 'Acetate', 321.0,
                          'Biomass_Ecol...', 0.8739215069684307],
                         ['D-Glucose', 10.0, 'H2O', 29.175827135565815, '',
                          ''],
                         ['Ammonium', 4.76531919319746, 'CO2',
                          22.80983331020498, '', ''],
                         ['Phosphate', 3.2148950476847533, 'H+',
                          17.530865429786687, '', '']])

    else:
        data = np.array([['o2_e', 21.799492655998794, 'ac_e', 321.0,
                          'Biomass_Ecol...', 0.8739215069684307],
                         ['glc__D_e', 10.0, 'h2o_e', 29.175827135565815, '',
                          ''],
                         ['nh4_e', 4.76531919319746, 'co2_e',
                          22.80983331020498, '', ''],
                         ['pi_e', 3.2148950476847533, 'h_e',
                          17.530865429786687, '', '']])

    expected_df = pd.DataFrame(
        data=data,
        columns=pd.MultiIndex.from_product(column_names)
    )

    assert expected_df.equals(out_df)


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
        data = np.array([['O2', 21.799492655998794, 'H2O', 29.175827135565815,
                          'Biomass_Ecol...', 0.8739215069684307],
                         ['D-Glucose', 10.0, 'CO2', 22.80983331020498, '',
                          ''],
                         ['Ammonium', 4.76531919319746, 'H+',
                          17.530865429786687, '', ''],
                         ['Phosphate', 3.2148950476847533, '', '', '', '']])

    else:
        data = np.array([['o2_e', 21.799492655998794, 'h2o_e',
                          29.175827135565815, 'Biomass_Ecol...',
                          0.8739215069684307],
                         ['glc__D_e', 10.0, 'co2_e', 22.80983331020498, '',
                          ''],
                         ['nh4_e', 4.76531919319746, 'h_e',
                          17.530865429786687, '', ''],
                         ['pi_e', 3.2148950476847533, '', '', '', '']])

    model.optimize()

    column_names = [['IN_FLUXES', 'OUT_FLUXES', 'OBJECTIVES'], ['ID', 'FLUX']]

    expected_df = pd.DataFrame(
        data=data,
        columns=pd.MultiIndex.from_product(column_names)
    )

    assert expected_df.equals(out_df)


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
    expected = np.array([['o2_e', 21.799493, 19.895962, 23.709518, 'h2o_e',
                          29.175827, 24.996702, 30.717036, 'Biomass_Ecol...',
                          0.8739215069684307],
                         ['glc__D_e', 10.0, 9.523306, 10.0, 'co2_e',
                          22.809833, 18.949008, 24.669342, '', ''],
                         ['nh4_e', 4.765319, 4.527053, 5.162646, 'h_e',
                          17.530865, 16.654322, 22.374655, '', ''],
                         ['pi_e', 3.214895, 3.05415, 3.214895, 'for_e', 0.0,
                          0.0, 5.720333, '', ''],
                         ['', '', '', '', 'ac_e', 0.0, 0.0, 1.906778, '', ''],
                         ['', '', '', '', 'acald_e', 0.0, 0.0, 1.271185, '',
                          ''],
                         ['', '', '', '', 'pyr_e', 0.0, 0.0, 1.271185, '',
                          ''],
                         ['', '', '', '', 'etoh_e', 0.0, 0.0, 1.107161, '',
                          ''],
                         ['', '', '', '', 'lac__D_e', 0.0, 0.0, 1.072563, '',
                          ''],
                         ['', '', '', '', 'succ_e', 0.0, 0.0, 0.837122, '',
                          ''],
                         ['', '', '', '', 'akg_e', 0.0, 0.0, 0.715042, '',
                          ''],
                         ['', '', '', '', 'glu__L_e', 0.0, 0.0, 0.635593, '',
                          '']])

    model.solver = opt_solver
    solution = model.optimize()
    out_df = model.summary(solution, fva=fraction).to_frame()

    assert np.array_equal(expected, out_df.values)


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

    expected = np.array([
        ['100%', 43.59898531199755,
         '2.0 h_c + 0.5 o2_c + q8h2_c --> h2o_c + 2.0 h_e...'],
        ['88%', 38.53460965051545,
         '4.0 h_c + nadh_c + q8_c --> 3.0 h_e + nad_c + q...'],
        ['12%', 5.064375661482105, 'q8_c + succ_c --> fum_c + q8h2_c']])

    out_df = model.metabolites.get_by_id(met).summary(solution).to_frame()

    assert np.array_equal(expected, out_df.values)


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

    out_df = model.metabolites.get_by_id(met).summary(names=names).to_frame()

    if names:
        data = np.array([
            ['100%', 43.59898531199755,
             '2.0 H+ + 0.5 O2 + Ubiquinol-8 --> H2O + 2.0 H+ ...'],
            ['88%', 38.53460965051545,
             '4.0 H+ + Nicotinamide adenine dinucleotide - re...'],
            ['12%', 5.064375661482105,
             'Ubiquinone-8 + Succinate --> Fumarate + Ubiquin...']])

    else:
        data = np.array([
            ['100%', 43.59898531199755,
             '2.0 h_c + 0.5 o2_c + q8h2_c --> h2o_c + 2.0 h_e...'],
            ['88%', 38.53460965051545,
             '4.0 h_c + nadh_c + q8_c --> 3.0 h_e + nad_c + q...'],
            ['12%', 5.064375661482105, 'q8_c + succ_c --> fum_c + q8h2_c']])

    assert np.array_equal(data, out_df.values)


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

    out_df = model.metabolites.get_by_id(met).summary(fva=fraction).to_frame()

    data = np.array([
        ['100%', 7.477382, 6.169728, 9.258708,
         'atp_c + f6p_c --> adp_c + fdp_c + h_c'],
        ['100%', 7.477382, 6.169728, 8.915488, 'fdp_c <=> dhap_c + g3p_c'],
        ['0%', 0.0, 0.0, 1.7161, 'fdp_c + h2o_c --> f6p_c + pi_c']])

    assert np.array(data, out_df.values)
