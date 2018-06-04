# -*- coding: utf-8 -*-

from __future__ import absolute_import, division

import logging
from operator import attrgetter

import pandas as pd
from numpy import zeros
from six import iteritems, print_
from six.moves import zip_longest
from tabulate import tabulate

from cobra.flux_analysis.variability import flux_variability_analysis
from cobra.util.solver import linear_reaction_coefficients
from cobra.util.util import format_long_string
from cobra.core import get_solution

LOGGER = logging.getLogger(__name__)


def metabolite_summary(met, solution=None, threshold=0.01, fva=False,
                       names=False, floatfmt='.3g'):
    """
    Print a summary of the production and consumption fluxes.

    This method requires the model for which this metabolite is a part
    to be solved.

    Parameters
    ----------
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
    if names:
        emit = attrgetter('name')
    else:
        emit = attrgetter('id')
    if solution is None:
        met.model.slim_optimize(error_value=None)
        solution = get_solution(met.model, reactions=met.reactions)

    rxns = sorted(met.reactions, key=attrgetter("id"))
    rxn_id = list()
    rxn_name = list()
    flux = list()
    reaction = list()
    for rxn in rxns:
        rxn_id.append(rxn.id)
        rxn_name.append(format_long_string(emit(rxn), 10))
        flux.append(solution[rxn.id] * rxn.metabolites[met])
        txt = rxn.build_reaction_string(use_metabolite_names=names)
        reaction.append(format_long_string(txt, 40 if fva is not None else 50))

    flux_summary = pd.DataFrame({
        "id": rxn_name,
        "flux": flux,
        "reaction": reaction
    }, index=rxn_id)

    if fva is not None:
        if hasattr(fva, 'columns'):
            fva_results = fva
        else:
            fva_results = flux_variability_analysis(
                met.model, list(met.reactions), fraction_of_optimum=fva)

        flux_summary["maximum"] = zeros(len(rxn_id), dtype=float)
        flux_summary["minimum"] = zeros(len(rxn_id), dtype=float)
        for rxn in rxns:
            fmax = rxn.metabolites[met] * fva_results.at[rxn.id, "maximum"]
            fmin = rxn.metabolites[met] * fva_results.at[rxn.id, "minimum"]
            if abs(fmin) <= abs(fmax):
                flux_summary.at[rxn.id, "fmax"] = fmax
                flux_summary.at[rxn.id, "fmin"] = fmin
            else:
                # Reverse fluxes.
                flux_summary.at[rxn.id, "fmax"] = fmin
                flux_summary.at[rxn.id, "fmin"] = fmax

    assert flux_summary["flux"].sum() < 1E-6, "Error in flux balance"

    flux_summary = _process_flux_dataframe(flux_summary, fva, threshold,
                                           floatfmt)

    flux_summary['percent'] = 0
    total_flux = flux_summary.loc[flux_summary.is_input, "flux"].sum()

    flux_summary.loc[flux_summary.is_input, 'percent'] = \
        flux_summary.loc[flux_summary.is_input, 'flux'] / total_flux
    flux_summary.loc[~flux_summary.is_input, 'percent'] = \
        flux_summary.loc[~flux_summary.is_input, 'flux'] / total_flux

    flux_summary['percent'] = flux_summary.percent.apply(
        lambda x: '{:.0%}'.format(x))

    if fva is not None:
        flux_table = tabulate(
            flux_summary.loc[:, ['percent', 'flux', 'fva_fmt', 'id',
                                 'reaction']].values, floatfmt=floatfmt,
            headers=['%', 'FLUX', 'RANGE', 'RXN ID', 'REACTION']).split('\n')
    else:
        flux_table = tabulate(
            flux_summary.loc[:, ['percent', 'flux', 'id', 'reaction']].values,
            floatfmt=floatfmt, headers=['%', 'FLUX', 'RXN ID', 'REACTION']
        ).split('\n')

    flux_table_head = flux_table[:2]

    met_tag = "{0} ({1})".format(format_long_string(met.name, 45),
                                 format_long_string(met.id, 10))

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


def model_summary(model, solution=None, threshold=0.01, fva=None, names=False,
                  floatfmt='.3g'):
    """
    Print a summary of the input and output fluxes of the model.

    Parameters
    ----------
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
    if names:
        emit = attrgetter('name')
    else:
        emit = attrgetter('id')
    objective_reactions = linear_reaction_coefficients(model)
    boundary_reactions = model.exchanges
    summary_rxns = set(objective_reactions.keys()).union(boundary_reactions)

    if solution is None:
        model.slim_optimize(error_value=None)
        solution = get_solution(model, reactions=summary_rxns)

    # Create a dataframe of objective fluxes
    obj_fluxes = pd.DataFrame({key: solution[key.id] * value for key,
                               value in iteritems(objective_reactions)},
                              index=['flux']).T
    obj_fluxes['id'] = obj_fluxes.apply(
        lambda x: format_long_string(x.name.id, 15), 1)

    # Build a dictionary of metabolite production from the boundary reactions
    metabolites = {m for r in boundary_reactions for m in r.metabolites}
    index = sorted(metabolites, key=attrgetter('id'))
    metabolite_fluxes = pd.DataFrame({
        'id': [format_long_string(emit(m), 15) for m in index],
        'flux': zeros(len(index), dtype=float)
    }, index=[m.id for m in index])
    for rxn in boundary_reactions:
        for met, stoich in iteritems(rxn.metabolites):
            metabolite_fluxes.at[met.id, 'flux'] += stoich * solution[rxn.id]

    # Calculate FVA results if requested
    if fva is not None:
        if len(index) != len(boundary_reactions):
            LOGGER.warning(
                "There exists more than one boundary reaction per metabolite. "
                "Please be careful when evaluating flux ranges.")
        metabolite_fluxes['fmin'] = zeros(len(index), dtype=float)
        metabolite_fluxes['fmax'] = zeros(len(index), dtype=float)
        if hasattr(fva, 'columns'):
            fva_results = fva
        else:
            fva_results = flux_variability_analysis(
                model, reaction_list=boundary_reactions,
                fraction_of_optimum=fva)

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
    metabolite_fluxes = _process_flux_dataframe(
        metabolite_fluxes, fva, threshold, floatfmt)

    # Begin building string output table
    def get_str_table(species_df, fva=False):
        """Formats a string table for each column"""
        if fva:
            return tabulate(
                species_df.loc[:, ['id', 'flux', 'fva_fmt']].values,
                floatfmt=floatfmt, tablefmt='simple',
                headers=['id', 'Flux', 'Range']).split('\n')
        else:
            return tabulate(species_df.loc[:, ['id', 'flux']].values,
                            floatfmt=floatfmt, tablefmt='plain').split('\n')

    in_table = get_str_table(
        metabolite_fluxes[metabolite_fluxes['is_input']], fva=fva is not None)
    out_table = get_str_table(
        metabolite_fluxes[~metabolite_fluxes['is_input']], fva=fva is not None)
    obj_table = get_str_table(obj_fluxes, fva=False)

    # Print nested output table
    print_(tabulate(
        [entries for entries in zip_longest(in_table, out_table, obj_table)],
        headers=['IN FLUXES', 'OUT FLUXES', 'OBJECTIVES'], tablefmt='simple'))


def _process_flux_dataframe(flux_dataframe, fva, threshold, floatfmt):
    """Some common methods for processing a database of flux information into
    print-ready formats. Used in both model_summary and metabolite_summary. """

    abs_flux = flux_dataframe['flux'].abs()
    flux_threshold = threshold * abs_flux.max()

    # Drop unused boundary fluxes
    if fva is None:
        flux_dataframe = flux_dataframe.loc[
            abs_flux >= flux_threshold, :].copy()
    else:
        flux_dataframe = flux_dataframe.loc[
            (abs_flux >= flux_threshold) |
            (flux_dataframe['fmin'].abs() >= flux_threshold) |
            (flux_dataframe['fmax'].abs() >= flux_threshold), :].copy()

        # Why set to zero? If included show true value?
        # flux_dataframe.loc[
        #     flux_dataframe['flux'].abs() < flux_threshold, 'flux'] = 0

    # Make all fluxes positive
    if fva is None:
        flux_dataframe['is_input'] = (flux_dataframe['flux'] >= 0)
        flux_dataframe['flux'] = flux_dataframe['flux'].abs()
    else:

        def get_direction(flux, fmin, fmax):
            """ decide whether or not to reverse a flux to make it positive """

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

    if fva is not None:
        flux_dataframe['fva_fmt'] = flux_dataframe.apply(
            lambda x: ("[{0.fmin:" + floatfmt + "}, {0.fmax:" +
                       floatfmt + "}]").format(x), 1)

        flux_dataframe = flux_dataframe.sort_values(
            by=['flux', 'fmax', 'fmin', 'id'],
            ascending=[False, False, False, True])

    else:
        flux_dataframe = flux_dataframe.sort_values(
            by=['flux', 'id'], ascending=[False, True])

    return flux_dataframe
