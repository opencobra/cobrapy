# -*- coding: utf-8 -*-

"""Define the ReactionSummary class."""

from __future__ import absolute_import

from operator import attrgetter

import pandas as pd
from six import iterkeys, itervalues

from cobra.core.summary import Summary


class ReactionSummary(Summary):
    """Define the ReactionSummary class.

    Parameters
    ----------
    rxn: cobra.Reaction
        The Reaction object whose summary we intend to get.
    names : bool
        Emit gene and metabolite names rather than identifiers.

    """

    def __init__(self, rxn, names, **kwargs):
        super(ReactionSummary, self).__init__(names=names, solution=None,
                                              threshold=None, fva=None,
                                              float_format=None, **kwargs)
        self.rxn = rxn

    def _generate(self):
        """
        Returns
        -------
        rxn_summary: pandas.DataFrame
            The DataFrame of reaction summary data.

        """
        if self.names:
            emit = attrgetter('name')
        else:
            emit = attrgetter('id')

        data = {
            'GENES_ID': [emit(gene) for gene in self.rxn.genes],
            'METABOLITES_ID':
            [emit(met) for met in iterkeys(self.rxn.metabolites)],
            'METABOLITES_STOICHIOMETRY':
            [met for met in itervalues(self.rxn.metabolites)],
            'METABOLITES_COMPARTMENT':
            [met.compartment for met in iterkeys(self.rxn.metabolites)]
        }

        rxn_summary = pd.DataFrame.from_dict(data, orient='index')\
                                  .T.fillna(value=pd.np.nan)

        rxn_summary.columns = pd.MultiIndex.from_tuples(
            [tuple(c.split('_')) for c in rxn_summary.columns])

        return rxn_summary

    def to_frame(self):
        """
        Returns
        -------
        A pandas.DataFrame of the summary.

        """
        return self._generate()

    def _to_table(self):
        """
        Returns
        -------
        A string of the summary table.

        """
        return self.to_frame().to_string(header=True, index=False, na_rep='',
                                         sparsify=False, justify='center')
