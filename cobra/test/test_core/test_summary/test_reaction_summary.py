# -*- coding: utf-8 -*-

"""Test functionalities of ReactionSummary."""

from __future__ import absolute_import

import numpy as np
import pytest

from cobra.test.test_core.test_summary import (
    captured_output, check_in_line, check_line)


@pytest.mark.skip(reason="has some weird formatting issue")
@pytest.mark.parametrize("rxn, names", [("ACALD", False), ("ACALD", True)])
def test_reaction_summary_to_table(model, rxn, names):
    """Test reaction summary._to_table()."""

    if names:
        expected_entry = [' adhE                                     '
                          'Coenzyme A             -1                   c    ']
    else:
        expected_entry = ['          accoa_c               '
                          '1                   c    ']

    with captured_output() as (out, _):
        print(model.reactions.get_by_id(rxn).summary(names=names))
    check_in_line(out.getvalue(), expected_entry)


@pytest.mark.parametrize("rxn, names", [("ACALD", False), ("FUM", True)])
def test_reaction_summary_to_frame(model, rxn, names):
    """Test reaction summary.to_frame()."""

    out_df = model.reactions.get_by_id(rxn).summary(names=names).to_frame()

    if names:
        expected_gene_names = ['fumB', 'fumC', 'fumA']
        expected_met_names = ['Fumarate', 'H2O', 'L-Malate']

    else:
        expected_gene_names = ['b0351', 'b1241', np.nan, np.nan, np.nan,
                               np.nan]
        expected_met_names = ['acald_c', 'coa_c', 'nad_c', 'accoa_c', 'h_c',
                              'nadh_c']

    assert all(out_df['GENES', 'ID'].tolist()) == \
        all(expected_gene_names)

    assert all(out_df['METABOLITES', 'ID'].tolist()) == \
        all(expected_met_names)
