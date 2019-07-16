# -*- coding: utf-8 -*-

"""Define the MetaboliteSummary class."""

from __future__ import absolute_import, division

from operator import attrgetter

import numpy as np
import pandas as pd

from cobra.core import get_solution
from cobra.core.summary import Summary
from cobra.flux_analysis.variability import flux_variability_analysis
from cobra.util import format_long_string


class MetaboliteSummary(Summary):
    """Define the MetaboliteSummary class.

    Parameters
    ----------
    met: cobra.Metabolite
        The Metabolite object whose summary we intend to get.
    solution : cobra.Solution
        A previously solved model solution to use for generating the
        summary. If None, the summary method will resolve the model.
        Note that the solution object must match the model, i.e., changes
        to the model such as changed bounds, added or removed reactions are
        not taken into account by this method.
    threshold : float
        Threshold below which fluxes are not reported.
    fva : pandas.DataFrame, float
        Whether or not to include flux variability analysis in the output.
        If given, fva should either be a previous FVA solution matching
        the model or a float between 0 and 1 representing the
        fraction of the optimum objective to be searched.
    names : bool
        Emit reaction and metabolite names rather than identifiers.
    float_format : one-parameter function
        Format string for floats.

    """

    def __init__(self, met, solution, threshold, fva, names, float_format,
                 **kwargs):
        super(MetaboliteSummary, self).__init__(solution, threshold, fva,
                                                names, float_format, **kwargs)
        self.met = met

    def _generate(self):
        """
        Returns
        -------
        flux_summary: pandas.DataFrame
            The DataFrame of flux summary data.

        """
        if self.names:
            emit = attrgetter('name')
        else:
            emit = attrgetter('id')

        if self.solution is None:
            self.met.model.slim_optimize(error_value=None)
            self.solution = get_solution(self.met.model,
                                         reactions=self.met.reactions)

        rxns = sorted(self.met.reactions, key=attrgetter('id'))

        data = [(emit(rxn), self.solution[rxn.id] * rxn.metabolites[self.met],
                 rxn.build_reaction_string(use_metabolite_names=self.names),
                 rxn) for rxn in rxns]

        flux_summary = pd.DataFrame.from_records(
            data=data,
            index=[rxn.id for rxn in rxns],
            columns=['id', 'flux', 'reaction_string', 'reaction']
        )

        assert flux_summary['flux'].sum() < self.met.model.tolerance, \
            'Error in flux balance'

        if self.fva is not None:
            if hasattr(self.fva, 'columns'):
                fva_results = self.fva
            else:
                fva_results = flux_variability_analysis(
                    self.met.model, list(self.met.reactions),
                    fraction_of_optimum=self.fva)

            flux_summary = pd.concat([flux_summary, fva_results],
                                     axis=1, sort=False)
            flux_summary.rename(columns={'maximum': 'fmax', 'minimum': 'fmin'},
                                inplace=True)

            def set_min_and_max(row):
                """Scale and set proper min and max values for flux."""
                fmax = row.reaction.metabolites[self.met] * row.fmax
                fmin = row.reaction.metabolites[self.met] * row.fmin

                if abs(fmin) <= abs(fmax):
                    row.fmax = fmax
                    row.fmin = fmin
                else:
                    # Reverse fluxes
                    row.fmax = fmin
                    row.fmin = fmax

                return row

            flux_summary = flux_summary.apply(set_min_and_max, axis=1)

        flux_summary = self._process_flux_dataframe(flux_summary)

        total_flux = flux_summary.loc[flux_summary.is_input, 'flux'].sum()

        # Calculate flux percentage
        flux_summary['percent'] = (flux_summary['flux'] / total_flux) * 100

        return flux_summary

    def to_frame(self):
        """
        Returns
        -------
        A pandas.DataFrame of the summary.

        """
        flux_df = self._generate()

        flux_df['is_input'] = flux_df['is_input']\
            .apply(lambda x: 'PRODUCING' if x is True else 'CONSUMING')

        flux_df.columns = [col.upper() for col in flux_df.columns]

        if self.fva is not None:
            flux_df.rename(columns={'IS_INPUT': 'RXN_STAT',
                                    'FMIN': 'FLUX_MIN',
                                    'FMAX': 'FLUX_MAX'}, inplace=True)
            flux_df = flux_df[['RXN_STAT', 'ID', 'PERCENT', 'FLUX', 'FLUX_MIN',
                               'FLUX_MAX', 'REACTION_STRING']]
        else:
            flux_df.rename(columns={'IS_INPUT': 'RXN_STAT'}, inplace=True)
            flux_df = flux_df[['RXN_STAT', 'ID', 'PERCENT', 'FLUX',
                               'REACTION_STRING']]

        flux_df.set_index(['RXN_STAT', 'ID'], inplace=True)

        return flux_df

    def _to_table(self):
        """
        Returns
        -------
        A string of the summary table.

        """
        return self.to_frame().to_string(header=True, index=True, na_rep='',
                                         float_format=self.float_format,
                                         sparsify=True, justify='center')
