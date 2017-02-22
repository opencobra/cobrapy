# -*- coding: utf-8 -*-
from __future__ import absolute_import

from warnings import warn

import numpy as np
from six import iteritems

try:
    import scipy
except ImportError:
    scipy = None


def create_stoichiometric_array(model, array_type='dense', dtype=None):
    """Return a stoichiometric array representation of the given model.

    The the columns represent the reactions and rows represent
    metabolites. S[i,j] therefore contains the quantity of metabolite `i`
    produced (negative for consumed) by reaction `j`.

    model: a :class:`~cobra.core.Model` object
    array_type: string
        The type of array to construct. if 'dense', return a standard
        numpy.array. Otherwise, 'dok', or 'lil' will construct a sparse array
        using scipy of the corresponding type.
    dtype: data-type
        The desired data-type for the array. If not given, defaults to float.

    """
    if array_type != 'dense' and not scipy:
        warn('Sparse matrices require scipy')

    if dtype is None:
        dtype = np.float64

    array_constructor = {
        'dense': np.zeros,
        'dok': scipy.sparse.dok_matrix,
        'lil': scipy.sparse.lil_matrix
    }

    n_metabolites = len(model.metabolites)
    n_reactions = len(model.reactions)
    array = array_constructor[array_type](
        (n_metabolites, n_reactions), dtype=dtype)

    # Convenience functions to index metabolites and reactions
    def m_ind(met): return model.metabolites.index(met)

    def r_ind(rxn): return model.reactions.index(rxn)

    for reaction in model.reactions:
        for metabolite, stoich in iteritems(reaction.metabolites):
            array[m_ind(metabolite), r_ind(reaction)] = stoich

    return array
