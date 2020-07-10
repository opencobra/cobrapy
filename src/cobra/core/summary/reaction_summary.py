# -*- coding: utf-8 -*-

"""Define the ReactionSummary class."""

from __future__ import absolute_import

from operator import attrgetter

import pandas as pd
from six import iterkeys, itervalues

from cobra.core.summary import Summary


class ReactionSummary(Summary):
    """
    Define the reaction summary.

    Attributes
    ----------
    reaction: cobra.Reaction
        The reaction to summarize.

    See Also
    --------
    Summary : Parent that defines further attributes.
    MetaboliteSummary
    ModelSummary

    """

    def __init__(self, reaction, model, **kwargs):
        """
        Initialize a metabolite summary.

        Parameters
        ----------
        reaction: cobra.Reaction
            The reaction object whose summary we intend to get.
        model : cobra.Model
            The metabolic model in which to generate a reaction summary.

        Other Parameters
        ----------------
        kwargs :
            Further keyword arguments are passed on to the parent class.

        See Also
        --------
        Summary : Parent that has further default parameters.
        MetaboliteSummary
        ModelSummary

        """
        super(ReactionSummary, self).__init__(model=model, **kwargs)
        self.reaction = reaction

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
            'GENES_ID': [emit(gene) for gene in self.reaction.genes],
            'METABOLITES_ID':
            [emit(met) for met in iterkeys(self.reaction.metabolites)],
            'METABOLITES_STOICHIOMETRY':
            [met for met in itervalues(self.reaction.metabolites)],
            'METABOLITES_COMPARTMENT':
            [met.compartment for met in iterkeys(self.reaction.metabolites)]
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
