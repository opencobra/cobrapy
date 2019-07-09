# -*- coding: utf-8 -*-

"""Define the ReactionSummary class."""

from __future__ import absolute_import

from operator import attrgetter

import pandas as pd
from six import iteritems

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

        gene_temp_df = pd.DataFrame([emit(gene) for gene in self.rxn.genes])
        met_temp_df = pd.DataFrame([
            [emit(key), value, key.compartment] for key, value in
            iteritems(self.rxn.metabolites)
        ])
        data = pd.concat([gene_temp_df, met_temp_df], axis=1).values

        columns = pd.MultiIndex.from_tuples((('REACTION', 'GENES', 'ID'),
                                             ('REACTION', 'METABOLITES', 'ID'),
                                             ('REACTION', 'METABOLITES',
                                              'STOICHIOMETRIC COEFFICIENT'),
                                             ('REACTION', 'METABOLITES',
                                              'COMPARTMENT')))

        rxn_summary = pd.DataFrame(
            data=data,
            columns=columns
        )

        del gene_temp_df, met_temp_df

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
