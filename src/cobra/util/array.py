# -*- coding: utf-8 -*-

from __future__ import absolute_import

from collections import namedtuple

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


def constraint_matrices(model, array_type='dense', include_vars=False,
                        zero_tol=1e-6):
    """Create a matrix representation of the problem.

    This is used for alternative solution approaches that do not use optlang.
    The function will construct the equality matrix, inequality matrix and
    bounds for the complete problem.

    Notes
    -----
    To accomodate non-zero equalities the problem will add the variable
    "const_one" which is a variable that equals one.

    Arguments
    ---------
    model : cobra.Model
        The model from which to obtain the LP problem.
    array_type : string
        The type of array to construct. if 'dense', return a standard
        numpy.array, 'dok', or 'lil' will construct a sparse array using
        scipy of the corresponding type and 'DataFrame' will give a
        pandas `DataFrame` with metabolite indices and reaction columns.
    zero_tol : float
        The zero tolerance used to judge whether two bounds are the same.

    Returns
    -------
    collections.namedtuple
        A named tuple consisting of 6 matrices and 2 vectors:
        - "equalities" is a matrix S such that S*vars = b. It includes a row
          for each constraint and one column for each variable.
        - "b" the right side of the equality equation such that S*vars = b.
        - "inequalities" is a matrix M such that lb <= M*vars <= ub.
          It contains a row for each inequality and as many columns as
          variables.
        - "bounds" is a compound matrix [lb ub] containing the lower and
          upper bounds for the inequality constraints in M.
        - "variable_fixed" is a boolean vector indicating whether the variable
          at that index is fixed (lower bound == upper_bound) and
          is thus bounded by an equality constraint.
        - "variable_bounds" is a compound matrix [lb ub] containing the
          lower and upper bounds for all variables.
    """
    if array_type not in ('DataFrame', 'dense') and not dok_matrix:
        raise ValueError('Sparse matrices require scipy')

    array_builder = {
        'dense': np.array, 'dok': dok_matrix, 'lil': lil_matrix,
        'DataFrame': pd.DataFrame,
    }[array_type]

    Problem = namedtuple("Problem",
                         ["equalities", "b", "inequalities", "bounds",
                          "variable_fixed", "variable_bounds"])
    equality_rows = []
    inequality_rows = []
    inequality_bounds = []
    b = []

    for const in model.constraints:
        lb = -np.inf if const.lb is None else const.lb
        ub = np.inf if const.ub is None else const.ub
        equality = (ub - lb) < zero_tol
        coefs = const.get_linear_coefficients(model.variables)
        coefs = [coefs[v] for v in model.variables]
        if equality:
            b.append(lb if abs(lb) > zero_tol else 0.0)
            equality_rows.append(coefs)
        else:
            inequality_rows.append(coefs)
            inequality_bounds.append([lb, ub])

    var_bounds = np.array([[v.lb, v.ub] for v in model.variables])
    fixed = var_bounds[:, 1] - var_bounds[:, 0] < zero_tol

    results = Problem(
        equalities=array_builder(equality_rows),
        b=np.array(b),
        inequalities=array_builder(inequality_rows),
        bounds=array_builder(inequality_bounds),
        variable_fixed=np.array(fixed),
        variable_bounds=array_builder(var_bounds))

    return results
