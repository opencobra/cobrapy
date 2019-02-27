# -*- coding: utf-8 -*-

"""Test functionalities of Geometric FBA."""

from __future__ import absolute_import

import numpy as np
import pytest
from pandas import Series

from cobra.core import Metabolite, Model, Reaction, Solution
from cobra.flux_analysis import geometric_fba


@pytest.fixture(scope="module")
def geometric_fba_model():
    """
    Generate geometric FBA model as described in [1]_

    References
    ----------
    .. [1] Smallbone, Kieran & Simeonidis, Vangelis. (2009).
           Flux balance analysis: A geometric perspective.
           Journal of theoretical biology.258. 311-5.
           10.1016/j.jtbi.2009.01.027.

    """
    test_model = Model('geometric_fba_paper_model')

    test_model.add_metabolites(Metabolite('A'))
    test_model.add_metabolites(Metabolite('B'))

    v_1 = Reaction('v1', upper_bound=1.0)
    v_1.add_metabolites({test_model.metabolites.A: 1.0})

    v_2 = Reaction('v2', lower_bound=-1000.0)
    v_2.add_metabolites({test_model.metabolites.A: -1.0,
                         test_model.metabolites.B: 1.0})

    v_3 = Reaction('v3', lower_bound=-1000.0)
    v_3.add_metabolites({test_model.metabolites.A: -1.0,
                         test_model.metabolites.B: 1.0})

    v_4 = Reaction('v4', lower_bound=-1000.0)
    v_4.add_metabolites({test_model.metabolites.A: -1.0,
                         test_model.metabolites.B: 1.0})

    v_5 = Reaction('v5')
    v_5.add_metabolites({test_model.metabolites.A: 0.0,
                         test_model.metabolites.B: -1.0})

    test_model.add_reactions([v_1, v_2, v_3, v_4, v_5])

    test_model.objective = 'v5'

    return test_model


def test_geometric_fba_benchmark(model, benchmark, all_solvers):
    """Benchmark geometric_fba."""
    model.solver = all_solvers
    benchmark(geometric_fba, model)


def test_geometric_fba(geometric_fba_model, all_solvers):
    """Test geometric_fba."""
    geometric_fba_model.solver = all_solvers
    geometric_fba_sol = geometric_fba(geometric_fba_model)
    expected = Series({'v1': 1.0, 'v2': 0.33, 'v3': 0.33, 'v4': 0.33,
                       'v5': 1.0}, index=['v1', 'v2', 'v3', 'v4', 'v5'])
    assert np.allclose(geometric_fba_sol.fluxes, expected, atol=1E-02)
