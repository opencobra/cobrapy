# -*- coding: utf-8 -*-
# cobra.flux_analysis.reaction.py
# functions for analyzing / creating objective functions

from __future__ import absolute_import

from six import iteritems

from cobra.core import Reaction


def assess(model, reaction, flux_coefficient_cutoff=0.001, solver=None):
    """Assesses the capacity of the model to produce the precursors for the
    reaction and absorb the production of the reaction while the reaction is
    operating at, or above, the specified cutoff.

    model: A :class:`~cobra.core.Model` object

    reaction: A :class:`~cobra.core.Reaction` object

    flux_coefficient_cutoff:  Float.  The minimum flux that reaction must carry
    to be considered active.

    solver : String or solver name. If None, the default solver will be used.

    returns: True if the model can produce the precursors and absorb the
    products for the reaction operating at, or above, flux_coefficient_cutoff.
    Otherwise, a dictionary of {'precursor': Status, 'product': Status}.  Where
    Status is the results from assess_precursors and assess_products,
    respectively.

    """
    reaction = model.reactions.get_by_id(reaction.id)
    original_objective = model.objective
    model.objective = reaction
    model.optimize(solver=solver)
    model.objective = original_objective
    if model.solution.f >= flux_coefficient_cutoff:
        return True
    else:
        results = {}
        results['precursors'] = assess_precursors(
            model, reaction, flux_coefficient_cutoff)
        results['products'] = assess_products(
            model, reaction, flux_coefficient_cutoff)
        return results


def assess_precursors(model, reaction, flux_coefficient_cutoff=0.001,
                      solver=None):
    """Assesses the ability of the model to provide sufficient precursors for
    a reaction operating at, or beyond, the specified cutoff.

    model: A :class:`~cobra.core.Model` object

    reaction: A :class:`~cobra.core.Reaction` object

    flux_coefficient_cutoff: Float. The minimum flux that reaction must carry
    to be considered active.

    solver : String or solver name. If None, the default solver will be used.

    returns: True if the precursors can be simultaneously produced at the
    specified cutoff. False, if the model has the capacity to produce each
    individual precursor at the specified threshold  but not all precursors at
    the required level simultaneously. Otherwise a dictionary of the required
    and the produced fluxes for each reactant that is not produced in
    sufficient quantities.

    """
    model = model.copy()
    reaction = model.reactions.get_by_id(reaction.id)
    original_objective = model.objective
    model.objective = reaction
    model.optimize(solver=solver)
    model.objective = original_objective
    if model.solution.f >= flux_coefficient_cutoff:
        return True
    #
    simulation_results = {}
    # build the sink reactions and add all at once
    sink_reactions = {}
    for the_component in reaction.reactants:
        # add in a sink reaction for each component
        sink_reaction = Reaction('test_sink_%s' % the_component.id)
        # then simulate production ability
        # then check it can exceed objective cutoff * component stoichiometric
        # coefficient.
        coefficient = reaction.get_coefficient(the_component)
        sink_reaction.add_metabolites({the_component: coefficient})
        sink_reaction.upper_bound = 1000
        sink_reactions[sink_reaction] = (the_component, coefficient)
    # First assess whether all precursors can pbe produced simultaneously
    super_sink = Reaction("super_sink")
    for reaction in sink_reactions:
        super_sink += reaction
    super_sink.id = 'super_sink'
    model.add_reactions(sink_reactions.keys() + [super_sink])
    model.objective = super_sink
    model.optimize(solver=solver)
    model.objective = original_objective
    if flux_coefficient_cutoff <= model.solution.f:
        return True

    # Otherwise assess the ability of the model to produce each precursor
    # individually.  Now assess the ability of the model to produce each
    # reactant for a reaction
    for sink_reaction, (component, coefficient) in iteritems(sink_reactions):
        # Calculate the maximum amount of the
        model.objective = sink_reaction
        model.optimize(solver=solver)
        model.objective = original_objective
        # metabolite that can be produced.
        if flux_coefficient_cutoff > model.solution.f:
            # Scale the results to a single unit
            simulation_results.update({
                component:
                    {
                        'required': flux_coefficient_cutoff / abs(coefficient),
                        'produced': model.solution.f / abs(coefficient)
                    }
            })
    if len(simulation_results) == 0:
        simulation_results = False
    return simulation_results


def assess_products(model, reaction, flux_coefficient_cutoff=0.001,
                    solver=None):
    """Assesses whether the model has the capacity to absorb the products of
    a reaction at a given flux rate.  Useful for identifying which components
    might be blocking a reaction from achieving a specific flux rate.

    model: A :class:`~cobra.core.Model` object

    reaction: A :class:`~cobra.core.Reaction` object

    flux_coefficient_cutoff:  Float.  The minimum flux that reaction must carry
    to be considered active.

    solver : String or solver name. If None, the default solver will be used.

    returns: True if the model has the capacity to absorb all the reaction
    products being simultaneously given the specified cutoff.   False, if the
    model has the capacity to absorb each individual product but not all
    products at the required level simultaneously.   Otherwise a dictionary of
    the required and the capacity fluxes for each product that is not absorbed
    in sufficient quantities.

    """
    model = model.copy()
    reaction = model.reactions.get_by_id(reaction.id)
    original_objective = model.objective
    model.objective = reaction
    model.optimize(solver=solver)
    model.objective = original_objective
    if model.solution.f >= flux_coefficient_cutoff:
        return True
    #
    simulation_results = {}
    # build the sink reactions and add all at once
    source_reactions = {}
    for the_component in reaction.products:
        # add in a sink reaction for each component
        source_reaction = Reaction('test_source_%s' % the_component.id)
        # then simulate production ability
        # then check it can exceed objective cutoff * component stoichiometric
        # coefficient.
        coefficient = reaction.get_coefficient(the_component)
        source_reaction.add_metabolites({the_component: coefficient})
        source_reaction.upper_bound = 1000
        source_reactions[source_reaction] = (the_component, coefficient)
    #
    super_source = Reaction('super_source')
    for reaction in source_reactions:
        super_source += reaction
    super_source.id = 'super_source'
    model.add_reactions(source_reactions.keys() + [super_source])
    model.objective = super_source
    model.optimize(solver=solver)
    model.objective = original_objective
    if flux_coefficient_cutoff <= model.solution.f:
        return True

    # Now assess the ability of the model to produce each reactant for a
    # reaction
    for source_reaction, (component, coefficient) in \
            iteritems(source_reactions):
        # Calculate the maximum amount of the
        model.objective = source_reaction
        model.optimize(solver=solver)
        model.objective = original_objective
        # metabolite that can be produced.
        if flux_coefficient_cutoff > model.solution.f:
            # Scale the results to a single unit
            simulation_results.update({
                component:
                    {
                        'required': flux_coefficient_cutoff / abs(coefficient),
                        'capacity': model.solution.f / abs(coefficient)}
                    }
            )
    if len(simulation_results) == 0:
        simulation_results = False
    return simulation_results
