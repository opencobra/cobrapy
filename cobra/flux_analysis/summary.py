# -*- coding: utf-8 -*-

from __future__ import absolute_import

import pandas as pd
from numpy import zeros
from six import iteritems, print_
from six.moves import zip_longest
from tabulate import tabulate

from cobra.flux_analysis.variability import flux_variability_analysis
from cobra.util.solver import linear_reaction_coefficients, choose_solver


def format_long_string(string, max_length):
    if len(string) > max_length:
        string = string[:max_length - 3]
        string += '...'
    return string


def metabolite_summary(met, threshold=0.01, fva=False, floatfmt='.3g',
                       **solver_args):
    """Print a summary of the reactions which produce and consume this
    metabolite

    threshold: float
    a value below which to ignore reaction fluxes

    fva: float (0->1), or None
    Whether or not to include flux variability analysis in the output.
    If given, fva should be a float between 0 and 1, representing the
    fraction of the optimum objective to be searched.

    floatfmt: string
        format method for floats, passed to tabulate. Default is '.3g'.

    """
    rxn_id = list()
    flux = list()
    reaction = list()
    for rxn in met.reactions:
        rxn_id.append(format_long_string(rxn.id, 10))
        flux.append(rxn.flux * rxn.metabolites[met])
        reaction.append(format_long_string(rxn.reaction, 40 if fva else 50))

    flux_summary = pd.DataFrame(data={
        "id": rxn_id, "flux": flux, "reaction": reaction})

    if fva:
        fva_results = flux_variability_analysis(
            met.model, met.reactions, fraction_of_optimum=fva,
            **solver_args)
        flux_summary.index = flux_summary["id"]
        flux_summary["maximum"] = zeros(len(rxn_id))
        flux_summary["minimum"] = zeros(len(rxn_id))
        for rid, rxn in zip(rxn_id, met.reactions):
            imax = rxn.metabolites[met] * fva_results.loc[rxn.id, "maximum"]
            imin = rxn.metabolites[met] * fva_results.loc[rxn.id, "minimum"]
            flux_summary.loc[rid, "fmax"] = (imax if abs(imin) <= abs(imax)
                                             else imin)
            flux_summary.loc[rid, "fmin"] = (imin if abs(imin) <= abs(imax)
                                             else imax)

    assert flux_summary.flux.sum() < 1E-6, "Error in flux balance"

    flux_summary = _process_flux_dataframe(flux_summary, fva, threshold,
                                           floatfmt)

    flux_summary['percent'] = 0
    total_flux = flux_summary[flux_summary.is_input].flux.sum()

    flux_summary.loc[flux_summary.is_input, 'percent'] = \
        flux_summary.loc[flux_summary.is_input, 'flux'] / total_flux
    flux_summary.loc[~flux_summary.is_input, 'percent'] = \
        flux_summary.loc[~flux_summary.is_input, 'flux'] / total_flux

    flux_summary['percent'] = flux_summary.percent.apply(
        lambda x: '{:.0%}'.format(x))

    if fva:
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


def model_summary(model, threshold=1E-8, fva=None, floatfmt='.3g',
                  **solver_args):
    """Print a summary of the input and output fluxes of the model.

    threshold: float
        tolerance for determining if a flux is zero (not printed)

    fva: int or None
        Whether or not to calculate and report flux variability in the
        output summary

    floatfmt: string
        format method for floats, passed to tabulate. Default is '.3g'.

    """
    legacy, _ = choose_solver(model, solver=solver_args.get("solver"))
    if legacy:
        raise NotImplementedError(
            "Summary support for legacy solvers was removed.")

    # Create a dataframe of objective fluxes
    objective_reactions = linear_reaction_coefficients(model)
    obj_fluxes = pd.DataFrame({key: key.flux * value for key, value in
                               iteritems(objective_reactions)},
                              index=['flux']).T
    obj_fluxes['id'] = obj_fluxes.apply(
        lambda x: format_long_string(x.name.id, 15), 1)

    # Build a dictionary of metabolite production from the boundary reactions
    # collect rxn.flux before fva which invalidates previous solver state
    boundary_reactions = model.exchanges
    metabolite_fluxes = {}
    for rxn in boundary_reactions:
        for met, stoich in iteritems(rxn.metabolites):
            metabolite_fluxes[met] = {
                'id': format_long_string(met.id, 15),
                'flux': stoich * rxn.flux}

    # Calculate FVA results if requested
    if fva:
        fva_results = flux_variability_analysis(
            model, reaction_list=boundary_reactions, fraction_of_optimum=fva,
            **solver_args)

        for rxn in boundary_reactions:
            for met, stoich in iteritems(rxn.metabolites):
                imin = stoich * fva_results.loc[rxn.id]['minimum']
                imax = stoich * fva_results.loc[rxn.id]['maximum']
                # Correct 'max' and 'min' for negative values
                metabolite_fluxes[met].update({
                    'fmin': imin if abs(imin) <= abs(imax) else imax,
                    'fmax': imax if abs(imin) <= abs(imax) else imin,
                })

    # Generate a dataframe of boundary fluxes
    metabolite_fluxes = pd.DataFrame(metabolite_fluxes).T
    metabolite_fluxes = _process_flux_dataframe(
        metabolite_fluxes, fva, threshold, floatfmt)

    # Begin building string output table
    def get_str_table(species_df, fva=False):
        """Formats a string table for each column"""

        if not fva:
            return tabulate(species_df.loc[:, ['id', 'flux']].values,
                            floatfmt=floatfmt, tablefmt='plain').split('\n')

        else:
            return tabulate(
                species_df.loc[:, ['id', 'flux', 'fva_fmt']].values,
                floatfmt=floatfmt, tablefmt='simple',
                headers=['id', 'Flux', 'Range']).split('\n')

    in_table = get_str_table(
        metabolite_fluxes[metabolite_fluxes.is_input], fva=fva)
    out_table = get_str_table(
        metabolite_fluxes[~metabolite_fluxes.is_input], fva=fva)
    obj_table = get_str_table(obj_fluxes, fva=False)

    # Print nested output table
    print_(tabulate(
        [entries for entries in zip_longest(in_table, out_table, obj_table)],
        headers=['IN FLUXES', 'OUT FLUXES', 'OBJECTIVES'], tablefmt='simple'))


def _process_flux_dataframe(flux_dataframe, fva, threshold, floatfmt):
    """Some common methods for processing a database of flux information into
    print-ready formats. Used in both model_summary and metabolite_summary. """

    # Drop unused boundary fluxes
    if not fva:
        flux_dataframe = flux_dataframe[
            flux_dataframe.flux.abs() > threshold].copy()
    else:
        flux_dataframe = flux_dataframe[
            (flux_dataframe.flux.abs() > threshold) |
            (flux_dataframe.fmin.abs() > threshold) |
            (flux_dataframe.fmax.abs() > threshold)].copy()

        flux_dataframe.loc[flux_dataframe.flux.abs() < threshold, 'flux'] = 0

    # Make all fluxes positive
    if not fva:
        flux_dataframe['is_input'] = flux_dataframe.flux >= 0
        flux_dataframe.flux = \
            flux_dataframe.flux.abs().astype('float')
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

    if fva:
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
