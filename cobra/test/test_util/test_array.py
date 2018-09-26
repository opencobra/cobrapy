# -*- coding: utf-8 -*-

"""Test functions of array.py"""

from __future__ import absolute_import

import numpy as np
import pytest

from cobra.util import create_stoichiometric_matrix

try:
    import scipy
except ImportError:
    scipy = None


def test_dense_matrix(model):
    S = create_stoichiometric_matrix(model, array_type='dense', dtype=int)
    assert S.dtype == int
    assert np.allclose(S.max(), [59])

    S_df = create_stoichiometric_matrix(
        model, array_type='DataFrame', dtype=int)
    assert S_df.values.dtype == int
    assert np.all(S_df.columns == [r.id for r in model.reactions])
    assert np.all(S_df.index == [m.id for m in model.metabolites])
    assert np.allclose(S_df.values, S)

    S = create_stoichiometric_matrix(model, array_type='dense', dtype=float)
    solution = model.optimize()
    mass_balance = S.dot(solution.fluxes)
    assert np.allclose(mass_balance, 0)


@pytest.mark.skipif(not scipy, reason='Sparse array methods require scipy')
def test_sparse_matrix(model):
    sparse_types = ['dok', 'lil']

    solution = model.optimize()
    for sparse_type in sparse_types:
        S = create_stoichiometric_matrix(model, array_type=sparse_type)
        mass_balance = S.dot(solution.fluxes)
        assert np.allclose(mass_balance, 0)
