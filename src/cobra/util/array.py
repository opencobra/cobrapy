"""Helper functions for array operations and sampling."""

from typing import TYPE_CHECKING, NamedTuple, Optional, Union

import numpy as np
import pandas as pd


# Used to avoid cyclic reference and enable third-party static type checkers to work
if TYPE_CHECKING:
    from cobra import Model


try:
    from scipy.sparse import dok_matrix, lil_matrix
except ImportError:
    dok_matrix, lil_matrix = None, None


def create_stoichiometric_matrix(
    model: "Model", array_type: str = "dense", dtype: Optional[np.dtype] = None
) -> Union[np.ndarray, dok_matrix, lil_matrix, pd.DataFrame]:
    """Return a stoichiometric array representation of the given model.

    The the columns represent the reactions and rows represent
    metabolites. S[i,j] therefore contains the quantity of metabolite `i`
    produced (negative for consumed) by reaction `j`.

    Parameters
    ----------
    model : cobra.Model
        The cobra model to construct the matrix for.
    array_type : {"dense", "dok", "lil", "DataFrame"}
        The type of array to construct. "dense" will return a standard
        numpy.ndarray. "dok", or "lil" will construct a sparse array using
        scipy of the corresponding type. "DataFrame" will give a
        pandas.DataFrame with metabolite as indices and reaction as
        columns.
    dtype : numpy.dtype, optional
        The desired numpy data type for the array (default numpy.float64).

    Returns
    -------
    matrix of class `dtype`
        The stoichiometric matrix for the given model.

    Raises
    ------
    ValueError
        If sparse matrix is used and scipy is not installed.

    .. deprecated:: 0.18.1
              "DataFrame" option for `array_type` will be replaced with
              "frame" in future versions.

    """
    if array_type not in ("DataFrame", "dense") and not dok_matrix:
        raise ValueError("Sparse matrices require scipy.")

    if dtype is None:
        dtype = np.float64

    array_constructor = {
        "dense": np.zeros,
        "dok": dok_matrix,
        "lil": lil_matrix,
        "DataFrame": np.zeros,
    }

    n_metabolites = len(model.metabolites)
    n_reactions = len(model.reactions)
    array = array_constructor[array_type]((n_metabolites, n_reactions), dtype=dtype)

    m_ind = model.metabolites.index
    r_ind = model.reactions.index

    for reaction in model.reactions:
        for metabolite, stoich in reaction.metabolites.items():
            array[m_ind(metabolite), r_ind(reaction)] = stoich

    if array_type == "DataFrame":
        metabolite_ids = [met.id for met in model.metabolites]
        reaction_ids = [rxn.id for rxn in model.reactions]
        return pd.DataFrame(array, index=metabolite_ids, columns=reaction_ids)

    else:
        return array


def nullspace(A: np.ndarray, atol: float = 1e-13, rtol: float = 0.0) -> np.ndarray:
    r"""Compute an approximate basis for the nullspace of A.

    The algorithm used by this function is based on the Singular Value
    Decomposition (SVD) of `A`.

    Parameters
    ----------
    A : numpy.ndarray
        `A` should be at most 2-D. 1-D array with length k will be treated
        as a 2-D with shape (1, k).
    atol : float, optional
        The absolute tolerance for a zero singular value. Singular values
        smaller than `atol` are considered to be zero (default 1e-13).
    rtol : float, optional
        The relative tolerance. Singular values less than `rtol * smax` are
        considered to be zero, where `smax` is the largest singular value
        (default 0.0).

    Returns
    -------
    numpy.ndarray
        If `A` is an array with shape (m, k), then `ns` will be an array
        with shape (k, n), where `n` is the estimated dimension of the
        nullspace of `A`.  The columns of `ns` are a basis for the
        nullspace; each element in numpy.dot(A, ns) will be approximately
        zero.

    Notes
    -----
    This is taken from the numpy cookbook.

    If both `atol` and `rtol` are positive, the combined tolerance is the
    maximum of the two; that is:

    .. math:: \mathtt{tol} = \max(\mathtt{atol}, \mathtt{rtol} * \mathtt{smax})

    Singular values smaller than `tol` are considered to be zero.

    """
    A = np.atleast_2d(A)
    _, s, vh = np.linalg.svd(A)
    tol = max(atol, rtol * s[0])
    nnz = (s >= tol).sum()
    ns = vh[nnz:].conj().T
    return ns


def constraint_matrices(
    model: "Model",
    array_type: str = "dense",
    zero_tol: float = 1e-6,
) -> NamedTuple:
    """Create a matrix representation of the problem.

    This is used for alternative solution approaches that do not use
    "optlang". The function will construct the equality matrix,
    inequality matrix and bounds for the complete problem.

    Parameters
    ----------
    model : cobra.Model
        The model from which to obtain the LP problem.
    array_type : {"dense", "dok", "lil", "DataFrame"}
        The type of array to construct. "dense" will return a standard
        numpy.ndarray. "dok", or "lil" will construct a sparse array using
        scipy of the corresponding type. "DataFrame" will give a
        pandas.DataFrame with metabolite as indices and reaction as
        columns.
    zero_tol : float, optional
        The zero tolerance used to judge whether two bounds are the same
        (default 1e-6).

    Returns
    -------
    NamedTuple
        A named tuple consisting of 6 matrices and 2 vectors:
        - "equalities" is a matrix `S` such that `S * vars = b`. It
          includes a row for each constraint and one column for each
          variable.
        - "b" is the right side of the equality equation such that
          `S * vars = b`.
        - "inequalities" is a matrix M such that `lb <= M * vars <= ub`.
          It contains a row for each inequality and as many columns as
          variables.
        - "bounds" is a compound matrix [lb ub] containing the lower and
          upper bounds for the inequality constraints in M.
        - "variable_fixed" is a boolean vector indicating whether the
          variable at that index is fixed (`lower bound == upper_bound`)
          and is thus bounded by an equality constraint.
        - "variable_bounds" is a compound matrix `[lb ub]` containing the
          lower and upper bounds for all variables.

    Notes
    -----
    To accomodate non-zero equalities, the problem will add the variable
    "const_one" which is a variable that equals one.

    .. deprecated:: 0.18.1
              "DataFrame" option for `array_type` will be replaced with
              "frame" in future versions.

    """
    if array_type not in ("DataFrame", "dense") and not dok_matrix:
        raise ValueError("Sparse matrices require scipy.")

    array_builder = {
        "dense": np.array,
        "dok": dok_matrix,
        "lil": lil_matrix,
        "DataFrame": pd.DataFrame,
    }[array_type]

    Problem = NamedTuple(
        "Problem",
        [
            ("equalities", Union[np.ndarray, dok_matrix, lil_matrix, pd.DataFrame]),
            ("b", np.ndarray),
            ("inequalities", Union[np.ndarray, dok_matrix, lil_matrix, pd.DataFrame]),
            ("bounds", Union[np.ndarray, dok_matrix, lil_matrix, pd.DataFrame]),
            ("variable_fixed", np.ndarray),
            (
                "variable_bounds",
                Union[np.ndarray, dok_matrix, lil_matrix, pd.DataFrame],
            ),
        ],
    )
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
        variable_bounds=array_builder(var_bounds),
    )

    return results
