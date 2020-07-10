# -*- coding: utf-8 -*-

"""Define the ModelSummary class."""

from __future__ import absolute_import

import logging
from collections import OrderedDict
from operator import attrgetter

import numpy as np
import pandas as pd
from six import iteritems, iterkeys

from cobra.core import get_solution
from cobra.core.summary import Summary
from cobra.flux_analysis.variability import flux_variability_analysis
from cobra.util.solver import linear_reaction_coefficients


logger = logging.getLogger(__name__)


class ModelSummary(Summary):
    """
    Define the model summary.

    See Also
    --------
    Summary : Parent that defines further attributes.
    MetaboliteSummary
    ReactionSummary

    """

    def __init__(self, model, **kwargs):
        """
        Initialize a model summary.

        Parameters
        ----------
        model : cobra.Model
            The metabolic model for which to generate a summary.

        Other Parameters
        ----------------
        kwargs :
            Further keyword arguments are passed on to the parent class.

        See Also
        --------
        Summary : Parent that has further default parameters.
        MetaboliteSummary
        ReactionSummary

        """
        super(ModelSummary, self).__init__(model=model, **kwargs)

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

        obj_rxns = linear_reaction_coefficients(self.model)
        boundary_rxns = self.model.exchanges
        summary_rxns = set(obj_rxns.keys()).union(boundary_rxns)

        if self.solution is None:
            self.model.slim_optimize(error_value=None)
            self.solution = get_solution(self.model, reactions=summary_rxns)

        # the order needs to be maintained
        ord_obj_rxns = OrderedDict(obj_rxns)

        obj_data = [(emit(rxn), self.solution[rxn.id] * stoich, np.nan, rxn,
                     np.nan, 1.0) for rxn, stoich in iteritems(ord_obj_rxns)]

        # create a DataFrame of objective reactions
        obj_df = pd.DataFrame.from_records(
            data=obj_data,
            index=[rxn.id for rxn in iterkeys(ord_obj_rxns)],
            columns=['id', 'flux', 'metabolite', 'reaction', 'is_input',
                     'is_obj_rxn']
        )

        # create a collection of metabolites from the boundary reactions
        mets = OrderedDict((met, rxn) for rxn in boundary_rxns
                           for met in sorted(rxn.metabolites,
                                             key=attrgetter('id')))
        index = [met.id for met in iterkeys(mets)]

        met_data = [(emit(met), self.solution[rxn.id] * rxn.metabolites[met],
                     met, rxn) for met, rxn in iteritems(mets)]

        # create a DataFrame of metabolites
        flux_summary = pd.DataFrame.from_records(
            data=met_data,
            index=index,
            columns=['id', 'flux', 'metabolite', 'reaction']
        )

        # Calculate FVA results if requested
        if self.fva is not None:
            if len(index) != len(boundary_rxns):
                logger.warning(
                    'There exists more than one boundary reaction per '
                    'metabolite. Please be careful when evaluating flux '
                    'ranges.')

            if hasattr(self.fva, 'columns'):
                fva_results = self.fva
            else:
                fva_results = flux_variability_analysis(
                    self.model, reaction_list=boundary_rxns,
                    fraction_of_optimum=self.fva
                )

            # save old index
            old_index = flux_summary.index
            # generate new index for consistency with fva_results
            flux_summary['rxn_id'] = flux_summary.apply(
                lambda x: x.reaction.id, axis=1
            )
            # change to new index
            flux_summary.set_index(['rxn_id'], inplace=True)
            flux_summary = pd.concat([flux_summary, fva_results],
                                     axis=1, sort=False)
            flux_summary.rename(columns={'maximum': 'fmax', 'minimum': 'fmin'},
                                inplace=True)
            # revert to old index
            flux_summary.set_index(old_index)

            def set_min_and_max(row):
                """Scale and set proper min and max values for flux."""
                fmax = row.reaction.metabolites[row.metabolite] * row.fmax
                fmin = row.reaction.metabolites[row.metabolite] * row.fmin

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

        flux_summary = pd.concat([flux_summary, obj_df], sort=False)

        return flux_summary

    def to_frame(self):
        """
        Returns
        -------
        A pandas.DataFrame of the summary.

        """
        flux_df = self._generate()

        # obtain separate DataFrames and join them instead of reorganizing the
        # old DataFrame, since that is easier to understand

        if self.fva is not None:
            column_names = ['IN_FLUXES-ID', 'IN_FLUXES-FLUX',
                            'IN_FLUXES-FLUX_MIN', 'IN_FLUXES-FLUX_MAX',
                            'OUT_FLUXES-ID', 'OUT_FLUXES-FLUX',
                            'OUT_FLUXES-FLUX_MIN', 'OUT_FLUXES-FLUX_MAX',
                            'OBJECTIVES-ID', 'OBJECTIVES-FLUX']

            # generate separate DataFrames
            met_in_df = flux_df[flux_df.is_input == 1.0]\
                .loc[:, ['id', 'flux', 'fmin', 'fmax']].reset_index(drop=True)

            met_out_df = flux_df[flux_df.is_input == 0.0]\
                .loc[:, ['id', 'flux', 'fmin', 'fmax']].reset_index(drop=True)

            obj_df = flux_df[flux_df.is_obj_rxn == 1.0]\
                .loc[:, ['id', 'flux']].reset_index(drop=True)

            # concatenate DataFrames
            concat_df = pd.concat([met_in_df, met_out_df, obj_df], axis=1)

            del met_in_df, met_out_df, obj_df

            # generate column names
            concat_df.columns = pd.MultiIndex.from_tuples(
                [tuple(c.split('-')) for c in column_names]
            )

        else:
            column_names = ['IN_FLUXES-ID', 'IN_FLUXES-FLUX',
                            'OUT_FLUXES-ID', 'OUT_FLUXES-FLUX',
                            'OBJECTIVES-ID', 'OBJECTIVES-FLUX']

            # generate separate DataFrames
            met_in_df = flux_df[flux_df.is_input == 1.0]\
                .loc[:, ['id', 'flux']].reset_index(drop=True)

            met_out_df = flux_df[flux_df.is_input == 0.0]\
                .loc[:, ['id', 'flux']].reset_index(drop=True)

            obj_df = flux_df[flux_df.is_obj_rxn == 1.0]\
                .loc[:, ['id', 'flux']].reset_index(drop=True)

            # concatenate DataFrames
            concat_df = pd.concat([met_in_df, met_out_df, obj_df], axis=1)

            del met_in_df, met_out_df, obj_df

            # generate column names
            concat_df.columns = pd.MultiIndex.from_tuples(
                [tuple(c.split('-')) for c in column_names]
            )

        return concat_df

    def _to_table(self):
        """
        Returns
        -------
        A string of the summary table.

        """
        return self.to_frame().to_string(header=True, index=False, na_rep='',
                                         float_format=self.float_format,
                                         sparsify=False, justify='center')
