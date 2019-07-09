# -*- coding: utf-8 -*-

"""Define the ModelSummary class."""

from __future__ import absolute_import

import logging
from operator import attrgetter

import numpy as np
import pandas as pd
from six import iteritems

from cobra.core import get_solution
from cobra.core.summary import Summary
from cobra.flux_analysis.variability import flux_variability_analysis
from cobra.util.solver import linear_reaction_coefficients
from cobra.util.util import format_long_string


logger = logging.getLogger(__name__)


class ModelSummary(Summary):
    """Define the ModelSummary class.

    Parameters
    ----------
    model: cobra.Model
        The Model object whose summary we intend to get.
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

    def __init__(self, model, solution, threshold, fva, names, float_format,
                 **kwargs):
        super(ModelSummary, self).__init__(solution, threshold, fva, names,
                                           float_format, **kwargs)
        self.model = model

    def _generate(self):
        """
        Returns
        -------
        metabolite_fluxes: pandas.DataFrame
            The DataFrame of metabolite fluxes.
        obj_fluxes: pandas.DataFrame
            The DataFrame of objective fluxes.

        """
        if self.names:
            emit = attrgetter('name')
        else:
            emit = attrgetter('id')

        objective_reactions = linear_reaction_coefficients(self.model)
        boundary_reactions = self.model.exchanges
        summary_rxns = \
            set(objective_reactions.keys()).union(boundary_reactions)

        if self.solution is None:
            self.model.slim_optimize(error_value=None)
            self.solution = get_solution(self.model, reactions=summary_rxns)

        # Create a pandas.DataFrame of objective fluxes
        obj_fluxes = pd.DataFrame({key: self.solution[key.id] * value
                                   for key, value in
                                   iteritems(objective_reactions)},
                                  index=['flux']).T
        obj_fluxes['id'] = obj_fluxes\
            .apply(lambda x: format_long_string(x.name.id, 15), axis=1)

        # Build a dictionary of metabolite production from the
        # boundary reactions
        metabolites = {met for rxn in boundary_reactions
                       for met in rxn.metabolites}
        index = sorted(metabolites, key=attrgetter('id'))

        # Create a pandas.DataFrame of metabolite fluxes
        metabolite_fluxes = pd.DataFrame({
            'id': [format_long_string(emit(met), 15) for met in index],
            'flux': np.zeros(len(index), dtype=float)
        }, index=[met.id for met in index])

        # Set the proper flux values for metabolites
        for rxn in boundary_reactions:
            for met, stoich in iteritems(rxn.metabolites):
                metabolite_fluxes.at[met.id, 'flux'] += \
                    stoich * self.solution[rxn.id]

        # Calculate FVA results if requested
        if self.fva is not None:
            if len(index) != len(boundary_reactions):
                logger.warning(
                    "There exists more than one boundary reaction per "
                    "metabolite. Please be careful when evaluating flux "
                    "ranges.")
            metabolite_fluxes['fmin'] = np.zeros(len(index), dtype=float)
            metabolite_fluxes['fmax'] = np.zeros(len(index), dtype=float)

            if hasattr(self.fva, 'columns'):
                fva_results = self.fva
            else:
                fva_results = flux_variability_analysis(
                    self.model, reaction_list=boundary_reactions,
                    fraction_of_optimum=self.fva)

            for rxn in boundary_reactions:
                for met, stoich in iteritems(rxn.metabolites):
                    fmin = stoich * fva_results.at[rxn.id, 'minimum']
                    fmax = stoich * fva_results.at[rxn.id, 'maximum']
                    # Correct 'max' and 'min' for negative values
                    if abs(fmin) <= abs(fmax):
                        metabolite_fluxes.at[met.id, 'fmin'] += fmin
                        metabolite_fluxes.at[met.id, 'fmax'] += fmax
                    else:
                        metabolite_fluxes.at[met.id, 'fmin'] += fmax
                        metabolite_fluxes.at[met.id, 'fmax'] += fmin

        # Generate a pandas.DataFrame of boundary fluxes
        metabolite_fluxes = self._process_flux_dataframe(metabolite_fluxes)

        return metabolite_fluxes, obj_fluxes

    def to_frame(self):
        """
        Returns
        -------
        A pandas.DataFrame of the summary.

        """
        met_df, obj_df = self._generate()

        if self.fva is not None:
            column_names = [['IN_FLUXES', 'OUT_FLUXES', 'OBJECTIVES'],
                            ['ID', 'FLUX', 'FLUX_MIN', 'FLUX_MAX']]

            # Generate new DataFrames for easy concatenation
            met_in_df = met_df[met_df['is_input']]\
                .loc[:, ['id', 'flux', 'fmin', 'fmax']].reset_index(drop=True)

            met_out_df = met_df[~met_df['is_input']]\
                .loc[:, ['id', 'flux', 'fmin', 'fmax']].reset_index(drop=True)

            obj_df = obj_df.loc[:, ['id', 'flux']].reset_index(drop=True)

            # concatenate DataFrames
            concat_df = pd.concat([met_in_df, met_out_df, obj_df], axis=1)\
                          .values

            del met_in_df, met_out_df, obj_df

            # Generate column names
            columns = pd.MultiIndex.from_product(column_names)

            # Remove 'FLUX_MIN' and 'FLUX_MAX' for 'OBJECTIVES'
            columns.set_codes([[0, 0, 0, 0, 2, 2, 2, 2, 1, 1],
                               [3, 0, 2, 1, 3, 0, 2, 1, 3, 0]],
                              inplace=True,
                              verify_integrity=False)

        else:
            column_names = [['IN_FLUXES', 'OUT_FLUXES', 'OBJECTIVES'],
                            ['ID', 'FLUX']]

            # Generate new DataFrames for easy concatenation
            met_in_df = met_df[met_df['is_input']]\
                .loc[:, ['id', 'flux']].reset_index(drop=True)

            met_out_df = met_df[~met_df['is_input']]\
                .loc[:, ['id', 'flux']].reset_index(drop=True)

            obj_df = obj_df.loc[:, ['id', 'flux']].reset_index(drop=True)

            # concatenate DataFrames
            concat_df = pd.concat([met_in_df, met_out_df, obj_df], axis=1)\
                          .values

            del met_in_df, met_out_df, obj_df

            # Generate column names
            columns = pd.MultiIndex.from_product(column_names)

        return pd.DataFrame(
            data=concat_df,
            columns=columns
        )

    def _to_table(self):
        """
        Returns
        -------
        A string of the summary table.

        """
        return self.to_frame().to_string(header=True, index=False, na_rep='',
                                         float_format=self.float_format,
                                         sparsify=False, justify='center')
