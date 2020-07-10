# -*- coding: utf-8 -*-

""" functions for analyzing / creating objective functions """

from __future__ import absolute_import, division

from operator import attrgetter
from warnings import warn

from six import iteritems

from cobra.core import Reaction


def assess(model, reaction, flux_coefficient_cutoff=0.001, solver=None):
    """Assesses production capacity.

    Assesses the capacity of the model to produce the precursors for the
    reaction and absorb the production of the reaction while the reaction is
    operating at, or above, the specified cutoff.

    Parameters
    ----------
    model : cobra.Model
        The cobra model to assess production capacity for

    reaction : reaction identifier or cobra.Reaction
        The reaction to assess

    flux_coefficient_cutoff :  float
        The minimum flux that reaction must carry to be considered active.

    solver : basestring
        Solver name. If None, the default solver will be used.

    Returns
    -------
    bool or dict
        True if the model can produce the precursors and absorb the products
        for the reaction operating at, or above, flux_coefficient_cutoff.
        Otherwise, a dictionary of {'precursor': Status, 'product': Status}.
        Where Status is the results from assess_precursors and
        assess_products, respectively.

    """
    reaction = model.reactions.get_by_any(reaction)[0]
    with model as m:
        m.objective = reaction
        if _optimize_or_value(m, solver=solver) >= flux_coefficient_cutoff:
            return True
        else:
            results = dict()
            results['precursors'] = assess_component(
                model, reaction, 'reactants', flux_coefficient_cutoff)
            results['products'] = assess_component(
                model, reaction, 'products', flux_coefficient_cutoff)
            return results


def assess_component(model, reaction, side, flux_coefficient_cutoff=0.001,
                     solver=None):
    """Assesses the ability of the model to provide sufficient precursors,
    or absorb products, for a reaction operating at, or beyond,
    the specified cutoff.

    Parameters
    ----------
    model : cobra.Model
        The cobra model to assess production capacity for

    reaction : reaction identifier or cobra.Reaction
        The reaction to assess

    side : basestring
        Side of the reaction, 'products' or 'reactants'

    flux_coefficient_cutoff :  float
        The minimum flux that reaction must carry to be considered active.

    solver : basestring
        Solver name. If None, the default solver will be used.

    Returns
    -------
    bool or dict
        True if the precursors can be simultaneously produced at the
        specified cutoff. False, if the model has the capacity to produce
        each individual precursor at the specified threshold  but not all
        precursors at the required level simultaneously. Otherwise a
        dictionary of the required and the produced fluxes for each reactant
        that is not produced in sufficient quantities.

    """
    reaction = model.reactions.get_by_any(reaction)[0]
    result_key = dict(reactants='produced', products='capacity')[side]
    get_components = attrgetter(side)
    with model as m:
        m.objective = reaction
        if _optimize_or_value(m, solver=solver) >= flux_coefficient_cutoff:
            return True
        simulation_results = {}
        # build the demand reactions and add all at once
        demand_reactions = {}
        for component in get_components(reaction):
            coeff = reaction.metabolites[component]
            demand = m.add_boundary(component, type='demand')
            demand.metabolites[component] = coeff
            demand_reactions[demand] = (component, coeff)
        # First assess whether all precursors can be produced simultaneously
        joint_demand = Reaction("joint_demand")
        for demand_reaction in demand_reactions:
            joint_demand += demand_reaction
        m.add_reactions([joint_demand])
        m.objective = joint_demand
        if _optimize_or_value(m, solver=solver) >= flux_coefficient_cutoff:
            return True

        # Otherwise assess the ability of the model to produce each precursor
        # individually.  Now assess the ability of the model to produce each
        # reactant for a reaction
        for demand_reaction, (component, coeff) in iteritems(demand_reactions):
            # Calculate the maximum amount of the
            with m:
                m.objective = demand_reaction
                flux = _optimize_or_value(m, solver=solver)
            # metabolite that can be produced.
            if flux_coefficient_cutoff > flux:
                # Scale the results to a single unit
                simulation_results.update({
                    component: {
                        'required': flux_coefficient_cutoff / abs(coeff),
                        result_key: flux / abs(coeff)
                    }})
        if len(simulation_results) == 0:
            simulation_results = False
        return simulation_results


def _optimize_or_value(model, value=0., solver=None):
    return model.slim_optimize(error_value=value)


def assess_precursors(model, reaction, flux_coefficient_cutoff=0.001,
                      solver=None):
    """Assesses the ability of the model to provide sufficient precursors for
    a reaction operating at, or beyond, the specified cutoff.

    Deprecated: use assess_component instead

    Parameters
    ----------
    model : cobra.Model
        The cobra model to assess production capacity for

    reaction : reaction identifier or cobra.Reaction
        The reaction to assess

    flux_coefficient_cutoff :  float
        The minimum flux that reaction must carry to be considered active.

    solver : basestring
        Solver name. If None, the default solver will be used.

    Returns
    -------
    bool or dict
        True if the precursors can be simultaneously produced at the
        specified cutoff. False, if the model has the capacity to produce
        each individual precursor at the specified threshold  but not all
        precursors at the required level simultaneously. Otherwise a
        dictionary of the required and the produced fluxes for each reactant
        that is not produced in sufficient quantities.

    """
    warn('use assess_component instead', DeprecationWarning)
    return assess_component(model, reaction, 'reactants',
                            flux_coefficient_cutoff, solver)


def assess_products(model, reaction, flux_coefficient_cutoff=0.001,
                    solver=None):
    """Assesses whether the model has the capacity to absorb the products of
    a reaction at a given flux rate.

    Useful for identifying which components might be blocking a reaction
    from achieving a specific flux rate.

    Deprecated: use assess_component instead

    Parameters
    ----------
    model : cobra.Model
        The cobra model to assess production capacity for

    reaction : reaction identifier or cobra.Reaction
        The reaction to assess

    flux_coefficient_cutoff :  float
        The minimum flux that reaction must carry to be considered active.

    solver : basestring
        Solver name. If None, the default solver will be used.

    Returns
    -------
    bool or dict
        True if the model has the capacity to absorb all the reaction
        products being simultaneously given the specified cutoff.   False,
        if the model has the capacity to absorb each individual product but
        not all products at the required level simultaneously.   Otherwise a
        dictionary of the required and the capacity fluxes for each product
        that is not absorbed in sufficient quantities.

    """
    warn('use assess_component instead', DeprecationWarning)
    return assess_component(model, reaction, 'products',
                            flux_coefficient_cutoff, solver)
