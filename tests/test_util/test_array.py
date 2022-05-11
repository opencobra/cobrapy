"""Test functions of array.py."""

from typing import TYPE_CHECKING

import numpy as np
import pytest

from cobra.util import create_stoichiometric_matrix


if TYPE_CHECKING:
    from cobra import Model


def test_dense_matrix(model: "Model") -> None:
    """Test dense stoichiometric matrix creation."""
    S = create_stoichiometric_matrix(model, array_type="dense", dtype=int)
    assert S.dtype == int
    assert np.allclose(S.max(), [59])

    S_df = create_stoichiometric_matrix(model, array_type="DataFrame", dtype=int)
    assert S_df.values.dtype == int
    assert np.all(S_df.columns == [r.id for r in model.reactions])
    assert np.all(S_df.index == [m.id for m in model.metabolites])
    assert np.allclose(S_df.values, S)

    S = create_stoichiometric_matrix(model, array_type="dense", dtype=float)
    solution = model.optimize()
    mass_balance = S.dot(solution.fluxes)
    assert np.allclose(mass_balance, 0)


def test_sparse_matrix(model: "Model") -> None:
    """Test sparse stoichiometric matrix creation."""
    pytest.importorskip("scipy")
    sparse_types = ["dok", "lil"]

    solution = model.optimize()
    for sparse_type in sparse_types:
        S = create_stoichiometric_matrix(model, array_type=sparse_type)
        mass_balance = S.dot(solution.fluxes)
        assert np.allclose(mass_balance, 0)
