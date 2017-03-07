# -*- coding: utf-8 -*-

from __future__ import absolute_import

import numpy as np
import pandas as pd
from six import iteritems

try:
    from scipy.sparse import dok_matrix, lil_matrix
except ImportError:
    dok_matrix, lil_matrix = None, None


def create_stoichiometric_matrix(model, array_type='dense', dtype=None):
    """Return a stoichiometric array representation of the given model.

    The the columns represent the reactions and rows represent
    metabolites. S[i,j] therefore contains the quantity of metabolite `i`
    produced (negative for consumed) by reaction `j`.

    Parameters
    ----------
    model : cobra.Model
        The cobra model to construct the matrix for.
    array_type : string
        The type of array to construct. if 'dense', return a standard
        numpy.array, 'dok', or 'lil' will construct a sparse array using
        scipy of the corresponding type and 'DataFrame' will give a
        pandas `DataFrame` with metabolite indices and reaction columns
    dtype : data-type
        The desired data-type for the array. If not given, defaults to float.

    Returns
    -------
    matrix of class `dtype`
        The stoichiometric matrix for the given model.
    """
    if array_type not in ('DataFrame', 'dense') and not dok_matrix:
        raise ValueError('Sparse matrices require scipy')

    if dtype is None:
        dtype = np.float64

    array_constructor = {
        'dense': np.zeros, 'dok': dok_matrix, 'lil': lil_matrix,
        'DataFrame': np.zeros,
    }

    n_metabolites = len(model.metabolites)
    n_reactions = len(model.reactions)
    array = array_constructor[array_type]((n_metabolites, n_reactions),
                                          dtype=dtype)

    m_ind = model.metabolites.index
    r_ind = model.reactions.index

    for reaction in model.reactions:
        for metabolite, stoich in iteritems(reaction.metabolites):
            array[m_ind(metabolite), r_ind(reaction)] = stoich

    if array_type == 'DataFrame':
        metabolite_ids = [met.id for met in model.metabolites]
        reaction_ids = [rxn.id for rxn in model.reactions]
        return pd.DataFrame(array, index=metabolite_ids, columns=reaction_ids)

    else:
        return array


def nullspace(A, atol=1e-13, rtol=0):
    """Compute an approximate basis for the nullspace of A.
    The algorithm used by this function is based on the singular value
    decomposition of `A`.

    Parameters
    ----------
    A : numpy.ndarray
        A should be at most 2-D.  A 1-D array with length k will be treated
        as a 2-D with shape (1, k)
    atol : float
        The absolute tolerance for a zero singular value.  Singular values
        smaller than `atol` are considered to be zero.
    rtol : float
        The relative tolerance.  Singular values less than rtol*smax are
        considered to be zero, where smax is the largest singular value.

    If both `atol` and `rtol` are positive, the combined tolerance is the
    maximum of the two; that is::
    tol = max(atol, rtol * smax)
    Singular values smaller than `tol` are considered to be zero.

    Returns
    -------
    numpy.ndarray
        If `A` is an array with shape (m, k), then `ns` will be an array
        with shape (k, n), where n is the estimated dimension of the
        nullspace of `A`.  The columns of `ns` are a basis for the
        nullspace; each element in numpy.dot(A, ns) will be approximately
        zero.

    Notes
    -----
    Taken from the numpy cookbook.
    """
    A = np.atleast_2d(A)
    u, s, vh = np.linalg.svd(A)
    tol = max(atol, rtol * s[0])
    nnz = (s >= tol).sum()
    ns = vh[nnz:].conj().T
    return ns
