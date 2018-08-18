# -*- coding: utf-8 -*-

"""Test functionalities of Geometric FBA."""

from __future__ import absolute_import

import numpy as np
from pandas import Series

from cobra.flux_analysis import geometric_fba


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
