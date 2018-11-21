# -*- coding: utf-8 -*-

from __future__ import absolute_import, division

from itertools import product

import pandas as pd
from numpy import abs, full, linspace, nan
from optlang.interface import OPTIMAL
from six import iteritems

import cobra.util.solver as sutil
from cobra.exceptions import OptimizationError
from cobra.flux_analysis import flux_variability_analysis as fva


def production_envelope(model, reactions, objective=None, carbon_sources=None,
                        points=20, threshold=1e-7):
    """Calculate the objective value conditioned on all combinations of
    fluxes for a set of chosen reactions

    The production envelope can be used to analyze a model's ability to
    produce a given compound conditional on the fluxes for another set of
    reactions, such as the uptake rates. The model is alternately optimized
    with respect to minimizing and maximizing the objective and the
    obtained fluxes are recorded. Ranges to compute production is set to the
    effective
    bounds, i.e., the minimum / maximum fluxes that can be obtained given
    current reaction bounds.

    Parameters
    ----------
    model : cobra.Model
        The model to compute the production envelope for.
    reactions : list or string
        A list of reactions, reaction identifiers or a single reaction.
    objective : string, dict, model.solver.interface.Objective, optional
        The objective (reaction) to use for the production envelope. Use the
        model's current objective if left missing.
    carbon_sources : list or string, optional
       One or more reactions or reaction identifiers that are the source of
       carbon for computing carbon (mol carbon in output over mol carbon in
       input) and mass yield (gram product over gram output). Only objectives
       with a carbon containing input and output metabolite is supported.
       Will identify active carbon sources in the medium if none are specified.
    points : int, optional
       The number of points to calculate production for.
    threshold : float, optional
        A cut-off under which flux values will be considered to be zero.

    Returns
    -------
    pandas.DataFrame
        A data frame with one row per evaluated point and

        - reaction id : one column per input reaction indicating the flux at
          each given point,
        - carbon_source: identifiers of carbon exchange reactions

        A column for the maximum and minimum each for the following types:

        - flux: the objective flux
        - carbon_yield: if carbon source is defined and the product is a
          single metabolite (mol carbon product per mol carbon feeding source)
        - mass_yield: if carbon source is defined and the product is a
          single metabolite (gram product per 1 g of feeding source)

    Examples
    --------
    >>> import cobra.test
    >>> from cobra.flux_analysis import production_envelope
    >>> model = cobra.test.create_test_model("textbook")
    >>> production_envelope(model, ["EX_glc__D_e", "EX_o2_e"])
    """
    reactions = model.reactions.get_by_any(reactions)
    objective = model.solver.objective if objective is None else objective
    data = dict()
    if carbon_sources is None:
        c_input = find_carbon_sources(model)
    else:
        c_input = model.reactions.get_by_any(carbon_sources)
    if c_input is None:
        data['carbon_source'] = None
    elif hasattr(c_input, 'id'):
        data['carbon_source'] = c_input.id
    else:
        data['carbon_source'] = ', '.join(rxn.id for rxn in c_input)

    size = points ** len(reactions)
    for direction in ('minimum', 'maximum'):
        data['flux_{}'.format(direction)] = full(size, nan, dtype=float)
        data['carbon_yield_{}'.format(direction)] = full(
            size, nan, dtype=float)
        data['mass_yield_{}'.format(direction)] = full(
            size, nan, dtype=float)
    grid = pd.DataFrame(data)
    with model:
        model.objective = objective
        objective_reactions = list(sutil.linear_reaction_coefficients(model))
        if len(objective_reactions) != 1:
            raise ValueError('cannot calculate yields for objectives with '
                             'multiple reactions')
        c_output = objective_reactions[0]
        min_max = fva(model, reactions, fraction_of_optimum=0)
        min_max[min_max.abs() < threshold] = 0.0
        points = list(product(*[
            linspace(min_max.at[rxn.id, "minimum"],
                     min_max.at[rxn.id, "maximum"],
                     points, endpoint=True) for rxn in reactions]))
        tmp = pd.DataFrame(points, columns=[rxn.id for rxn in reactions])
        grid = pd.concat([grid, tmp], axis=1, copy=False)
        add_envelope(model, reactions, grid, c_input, c_output, threshold)
    return grid


def add_envelope(model, reactions, grid, c_input, c_output, threshold):
    if c_input is not None:
        input_components = [reaction_elements(rxn) for rxn in c_input]
        output_components = reaction_elements(c_output)
        try:
            input_weights = [reaction_weight(rxn) for rxn in c_input]
            output_weight = reaction_weight(c_output)
        except ValueError:
            input_weights = []
            output_weight = []
    else:
        input_components = []
        output_components = []
        input_weights = []
        output_weight = []

    for direction in ('minimum', 'maximum'):
        with model:
            model.objective_direction = direction
            for i in range(len(grid)):
                with model:
                    for rxn in reactions:
                        point = grid.at[i, rxn.id]
                        rxn.bounds = point, point
                    obj_val = model.slim_optimize()
                    if model.solver.status != OPTIMAL:
                        continue
                    grid.at[i, 'flux_{}'.format(direction)] = \
                        0.0 if abs(obj_val) < threshold else obj_val
                    if c_input is not None:
                        grid.at[i, 'carbon_yield_{}'.format(direction)] = \
                            total_yield([rxn.flux for rxn in c_input],
                                        input_components,
                                        obj_val,
                                        output_components)
                        grid.at[i, 'mass_yield_{}'.format(direction)] = \
                            total_yield([rxn.flux for rxn in c_input],
                                        input_weights,
                                        obj_val,
                                        output_weight)


def total_yield(input_fluxes, input_elements, output_flux, output_elements):
    """
    Compute total output per input unit.

    Units are typically mol carbon atoms or gram of source and product.

    Parameters
    ----------
    input_fluxes : list
        A list of input reaction fluxes in the same order as the
        ``input_components``.
    input_elements : list
        A list of reaction components which are in turn list of numbers.
    output_flux : float
        The output flux value.
    output_elements : list
        A list of stoichiometrically weighted output reaction components.

    Returns
    -------
    float
        The ratio between output (mol carbon atoms or grams of product) and
        input (mol carbon atoms or grams of source compounds).
    """

    carbon_input_flux = sum(
        total_components_flux(flux, components, consumption=True)
        for flux, components in zip(input_fluxes, input_elements))
    carbon_output_flux = total_components_flux(
        output_flux, output_elements, consumption=False)
    try:
        return carbon_output_flux / carbon_input_flux
    except ZeroDivisionError:
        return nan


def reaction_elements(reaction):
    """
    Split metabolites into the atoms times their stoichiometric coefficients.

    Parameters
    ----------
    reaction : Reaction
        The metabolic reaction whose components are desired.

    Returns
    -------
    list
        Each of the reaction's metabolites' desired carbon elements (if any)
        times that metabolite's stoichiometric coefficient.
    """
    c_elements = [coeff * met.elements.get('C', 0)
                  for met, coeff in iteritems(reaction.metabolites)]
    return [elem for elem in c_elements if elem != 0]


def reaction_weight(reaction):
    """Return the metabolite weight times its stoichiometric coefficient."""
    if len(reaction.metabolites) != 1:
        raise ValueError('Reaction weight is only defined for single '
                         'metabolite products or educts.')
    met, coeff = next(iteritems(reaction.metabolites))
    return [coeff * met.formula_weight]


def total_components_flux(flux, components, consumption=True):
    """
    Compute the total components consumption or production flux.

    Parameters
    ----------
    flux : float
        The reaction flux for the components.
    components : list
        List of stoichiometrically weighted components.
    consumption : bool, optional
        Whether to sum up consumption or production fluxes.

    """
    direction = 1 if consumption else -1
    c_flux = [elem * flux * direction for elem in components]
    return sum([flux for flux in c_flux if flux > 0])


def find_carbon_sources(model):
    """
    Find all active carbon source reactions.

    Parameters
    ----------
    model : Model
        A genome-scale metabolic model.

    Returns
    -------
    list
       The medium reactions with carbon input flux.

    """
    try:
        model.slim_optimize(error_value=None)
    except OptimizationError:
        return []

    reactions = model.reactions.get_by_any(list(model.medium))
    reactions_fluxes = [
        (rxn, total_components_flux(rxn.flux, reaction_elements(rxn),
                                    consumption=True)) for rxn in reactions]
    return [rxn for rxn, c_flux in reactions_fluxes if c_flux > 0]
