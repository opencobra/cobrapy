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

        rxns = sorted(self.met.reactions, key=attrgetter("id"))
        rxn_id = []
        rxn_name = []
        flux = []
        reaction = []

        for rxn in rxns:
            rxn_id.append(rxn.id)
            rxn_name.append(format_long_string(emit(rxn), 10))
            flux.append(self.solution[rxn.id] * rxn.metabolites[self.met])
            txt = rxn.build_reaction_string(use_metabolite_names=self.names)
            reaction.append(format_long_string(txt, 40
                                               if self.fva is not None
                                               else 50))

        flux_summary = pd.DataFrame({
            "id": rxn_name,
            "flux": flux,
            "reaction": reaction
        }, index=rxn_id)

        if self.fva is not None:
            if hasattr(self.fva, 'columns'):
                fva_results = self.fva
            else:
                fva_results = flux_variability_analysis(
                    self.met.model, list(self.met.reactions),
                    fraction_of_optimum=self.fva)

            flux_summary["maximum"] = np.zeros(len(rxn_id), dtype=float)
            flux_summary["minimum"] = np.zeros(len(rxn_id), dtype=float)

            for rxn in rxns:
                fmax = rxn.metabolites[self.met] * \
                    fva_results.at[rxn.id, "maximum"]

                fmin = rxn.metabolites[self.met] * \
                    fva_results.at[rxn.id, "minimum"]

                if abs(fmin) <= abs(fmax):
                    flux_summary.at[rxn.id, "fmax"] = fmax
                    flux_summary.at[rxn.id, "fmin"] = fmin
                else:
                    # Reverse fluxes.
                    flux_summary.at[rxn.id, "fmax"] = fmin
                    flux_summary.at[rxn.id, "fmin"] = fmax

        assert flux_summary["flux"].sum() < 1E-6, "Error in flux balance"

        flux_summary = self._process_flux_dataframe(flux_summary)

        flux_summary['percent'] = 0
        total_flux = flux_summary.loc[flux_summary.is_input, "flux"].sum()

        flux_summary.loc[flux_summary.is_input, 'percent'] = \
            flux_summary.loc[flux_summary.is_input, 'flux'] / total_flux

        flux_summary.loc[~flux_summary.is_input, 'percent'] = \
            flux_summary.loc[~flux_summary.is_input, 'flux'] / total_flux

        flux_summary['percent'] = flux_summary.percent.\
            apply(lambda x: '{:.0%}'.format(x))

        return flux_summary

    def to_frame(self):
        """
        Returns
        -------
        A pandas.DataFrame of the summary.

        """
        flux_df = self._generate()

        if self.fva is not None:
            column_names = ['percent', 'flux', 'fmin',
                            'fmax', 'id', 'reaction']

            # Generate DataFrame of production reactions
            flux_prod = flux_df[flux_df.is_input.values]\
                .loc[:, column_names]\
                .reset_index(drop=True)
            flux_prod.columns = [name.upper() for name in column_names]
            flux_prod['RXN_STAT'] = 'PRODUCING'

            # Generate DataFrame of consumption reactions
            flux_cons = flux_df[~flux_df.is_input.values]\
                .loc[:, column_names]\
                .reset_index(drop=True)
            flux_cons.columns = [name.upper() for name in column_names]
            flux_cons['RXN_STAT'] = 'CONSUMING'

            concat_df = pd.concat([flux_prod, flux_cons])

            del flux_prod, flux_cons

            concat_df.rename(columns={'FMIN': 'FLUX_MIN', 'FMAX': 'FLUX_MAX'},
                             inplace=True)
            concat_df.set_index(['RXN_STAT', 'ID'], inplace=True)

        else:
            column_names = ['percent', 'flux', 'id', 'reaction']

            # Generate DataFrame of production reactions
            flux_prod = flux_df[flux_df.is_input.values]\
                .loc[:, column_names]\
                .reset_index(drop=True)
            flux_prod.columns = [name.upper() for name in column_names]
            flux_prod['RXN_STAT'] = 'PRODUCING'

            # Generate DataFrame of consumption reactions
            flux_cons = flux_df[~flux_df.is_input.values]\
                .loc[:, column_names]\
                .reset_index(drop=True)
            flux_cons.columns = [name.upper() for name in column_names]
            flux_cons['RXN_STAT'] = 'CONSUMING'

            concat_df = pd.concat([flux_prod, flux_cons])

            del flux_prod, flux_cons

            concat_df.set_index(['RXN_STAT', 'ID'], inplace=True)

        return concat_df

    def _to_table(self):
        """
        Returns
        -------
        A string of the summary table.

        """
        return self.to_frame().to_string(header=True, index=True, na_rep='',
                                         float_format=self.float_format,
                                         sparsify=True, justify='center')
