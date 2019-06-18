# -*- coding: utf-8 -*-

"""Define the Summary class."""

from __future__ import absolute_import, division

import logging
from operator import attrgetter

import pandas as pd
import numpy as np
from six import iteritems, print_
from six.moves import zip_longest
from tabulate import tabulate

from cobra.core import get_solution
from cobra.flux_analysis.variability import flux_variability_analysis
from cobra.util.solver import linear_reaction_coefficients
from cobra.util.util import format_long_string


LOGGER = logging.getLogger(__name__)


class Summary(object):
    """Class definition for a Summary object.

    Parameters
    ----------
    solution : cobra.Solution or None
        A previously solved model solution to use for generating the
        summary. If none provided (default), the summary method will
        resolve the model. Note that the solution object must match the
        model, i.e., changes to the model such as changed bounds,
        added or removed reactions are not taken into account by this
        method.
    threshold : float
        Threshold below which fluxes are not reported.
    fva : pandas.DataFrame, float or None
        Whether or not to include flux variability analysis in the output.
        If given, fva should either be a previous FVA solution matching
        the model or a float between 0 and 1 representing the
        fraction of the optimum objective to be searched.
    names : bool
        Emit reaction and metabolite names rather than identifiers (default
        False).
    floatfmt : string
        Format string for floats (default '.3g').

    """

    def __init__(self, solution, threshold, fva, names, floatfmt):
        self.solution = solution
        self.threshold = threshold
        self.fva = fva
        self.names = names
        self.floatfmt = floatfmt

    def generate():
        """Generate the summary for the required cobra object."""
        pass

    def _process_flux_dataframe(self, flux_dataframe):
        """Some common methods for processing a database of flux information
        into print-ready formats. Used in both ModelSummary and
        MetaboliteSummary."""

        abs_flux = flux_dataframe['flux'].abs()
        flux_threshold = self.threshold * abs_flux.max()

        # Drop unused boundary fluxes
        if self.fva is None:
            flux_dataframe = \
                flux_dataframe.loc[abs_flux >= flux_threshold, :].copy()
        else:
            flux_dataframe = (
                flux_dataframe
                .loc[(abs_flux >= flux_threshold) |
                     (flux_dataframe['fmin'].abs() >= flux_threshold) |
                     (flux_dataframe['fmax'].abs() >= flux_threshold), :]
                .copy()
                )

            # Why set to zero? If included show true value?
            # flux_dataframe.loc[
            #     flux_dataframe['flux'].abs() < flux_threshold, 'flux'] = 0

        # Make all fluxes positive
        if self.fva is None:
            flux_dataframe['is_input'] = (flux_dataframe['flux'] >= 0)
            flux_dataframe['flux'] = flux_dataframe['flux'].abs()
        else:

            def get_direction(flux, fmin, fmax):
                """Decide whether or not to reverse a flux to make it
                positive."""

                if flux < 0:
                    return -1
                elif flux > 0:
                    return 1
                elif (fmax > 0) & (fmin <= 0):
                    return 1
                elif (fmax < 0) & (fmin >= 0):
                    return -1
                elif ((fmax + fmin) / 2) < 0:
                    return -1
                else:
                    return 1

            sign = flux_dataframe.apply(
                lambda x: get_direction(x.flux, x.fmin, x.fmax), 1)

            flux_dataframe['is_input'] = sign == 1

            flux_dataframe.loc[:, ['flux', 'fmin', 'fmax']] = \
                flux_dataframe.loc[:, ['flux', 'fmin', 'fmax']].multiply(
                    sign, 0).astype('float').round(6)

            flux_dataframe.loc[:, ['flux', 'fmin', 'fmax']] = \
                flux_dataframe.loc[:, ['flux', 'fmin', 'fmax']].applymap(
                    lambda x: x if abs(x) > 1E-6 else 0)

        if self.fva is not None:
            flux_dataframe['fva_fmt'] = flux_dataframe.apply(
                lambda x: ("[{0.fmin:" + self.floatfmt + "},\
                {0.fmax:" + self.floatfmt + "}]").format(x), 1)

            flux_dataframe = flux_dataframe.sort_values(
                by=['flux', 'fmax', 'fmin', 'id'],
                ascending=[False, False, False, True])

        else:
            flux_dataframe = flux_dataframe.sort_values(
                by=['flux', 'id'], ascending=[False, True])

        return flux_dataframe

    def _repr_html_(self):
        return self.generate()


class MetaboliteSummary(Summary):
    """Class definition for a MetaboliteSummary object.

    Parameters
    ----------
    met: cobra.Metabolite
        The Metabolite object whose summary we intend to get.
    solution : cobra.Solution, optional
        A previously solved model solution to use for generating the
        summary. If none provided (default), the summary method will
        resolve the model. Note that the solution object must match the
        model, i.e., changes to the model such as changed bounds,
        added or removed reactions are not taken into account by this
        method.
    threshold : float, optional
        Threshold below which fluxes are not reported.
    fva : pandas.DataFrame, float or None, optional
        Whether or not to include flux variability analysis in the output.
        If given, fva should either be a previous FVA solution matching
        the model or a float between 0 and 1 representing the
        fraction of the optimum objective to be searched.
    names : bool, optional
        Emit reaction and metabolite names rather than identifiers (default
        False).
    floatfmt : string, optional
        Format string for floats (default '.3g').

    """

    def __init__(self, met, solution, threshold, fva, names, floatfmt):
        super(MetaboliteSummary, self).__init__(solution, threshold, fva,
                                                names, floatfmt)
        self.met = met

    def generate(self):
        if self.names:
            emit = attrgetter('name')
        else:
            emit = attrgetter('id')

        if self.solution is None:
            self.met.model.slim_optimize(error_value=None)
            self.solution = get_solution(self.met.model,
                                         reactions=self.met.reactions)

        rxns = sorted(self.met.reactions, key=attrgetter("id"))
        rxn_id = list()
        rxn_name = list()
        flux = list()
        reaction = list()

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

        flux_summary['percent'] = flux_summary.percent.apply(
            lambda x: '{:.0%}'.format(x))

        if self.fva is not None:
            flux_table = tabulate(
                flux_summary.loc[:, ['percent', 'flux', 'fva_fmt', 'id',
                                     'reaction']].values,
                floatfmt=self.floatfmt,
                headers=['%', 'FLUX', 'RANGE', 'RXN ID',
                         'REACTION']).split('\n')
        else:
            flux_table = tabulate(
                flux_summary.loc[:, ['percent', 'flux', 'id',
                                     'reaction']].values,
                floatfmt=self.floatfmt,
                headers=['%', 'FLUX', 'RXN ID', 'REACTION']
            ).split('\n')

        flux_table_head = flux_table[:2]

        met_tag = "{0} ({1})".format(format_long_string(self.met.name, 45),
                                     format_long_string(self.met.id, 10))

        head = "PRODUCING REACTIONS -- " + met_tag
        print_(head)
        print_("-" * len(head))
        print_('\n'.join(flux_table_head))
        print_('\n'.join(
            pd.np.array(flux_table[2:])[flux_summary.is_input.values]))

        print_()
        print_("CONSUMING REACTIONS -- " + met_tag)
        print_("-" * len(head))
        print_('\n'.join(flux_table_head))
        print_('\n'.join(
            pd.np.array(flux_table[2:])[~flux_summary.is_input.values]))


class ModelSummary(Summary):
    """Class definition for a ModelSummary object.

    Parameters
    ----------
    model: cobra.Model
        The Model object whose summary we intend to get.
    solution: cobra.Solution, optional
        A previously solved model solution to use for generating the
        summary. If none provided (default), the summary method will
        resolve the model. Note that the solution object must match the
        model, i.e., changes to the model such as changed bounds,
        added or removed reactions are not taken into account by this
        method.
    threshold : float, optional
        Threshold below which fluxes are not reported.
    fva : pandas.DataFrame, float or None, optional
        Whether or not to include flux variability analysis in the output.
        If given, fva should either be a previous FVA solution matching
        the model or a float between 0 and 1 representing the
        fraction of the optimum objective to be searched.
    names : bool, optional
        Emit reaction and metabolite names rather than identifiers (default
        False).
    floatfmt : string, optional
        Format string for floats (default '.3g').

    """

    def __init__(self, model, solution, threshold, fva, names, floatfmt):
        super(ModelSummary, self).__init__(solution, threshold, fva, names,
                                           floatfmt)
        self.model = model

    def generate(self):
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
                LOGGER.warning(
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

        # Generate a dataframe of boundary fluxes
        metabolite_fluxes = self._process_flux_dataframe(metabolite_fluxes)

        # Begin building string output table
        def get_str_table(species_df, fva=False):
            """Formats a string table for each column"""
            if fva:
                return tabulate(
                    species_df.loc[:, ['id', 'flux', 'fva_fmt']].values,
                    floatfmt=self.floatfmt, tablefmt='simple',
                    headers=['id', 'Flux', 'Range']).split('\n')
            else:
                return tabulate(species_df.loc[:, ['id', 'flux']].values,
                                floatfmt=self.floatfmt,
                                tablefmt='plain').split('\n')

        in_table = get_str_table(met_df[met_df['is_input']],
                                 fva=self.fva is not None)

        out_table = get_str_table(met_df[~met_df['is_input']],
                                  fva=self.fva is not None)

        obj_table = get_str_table(obj_df, fva=False)

        # Print nested output table
        print_(tabulate(
            [entries for entries in zip_longest(in_table, out_table,
                                                obj_table)],
            headers=['IN FLUXES', 'OUT FLUXES', 'OBJECTIVES'],
            tablefmt='simple'))
