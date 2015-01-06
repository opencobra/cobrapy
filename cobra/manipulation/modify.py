from copy import deepcopy
from warnings import warn

from .. import Reaction


def initialize_growth_medium(cobra_model, the_medium='MgM',
                             external_boundary_compartment='e',
                             external_boundary_reactions=None,
                             reaction_lower_bound=0.,
                             reaction_upper_bound=1000.,
                             irreversible=False,
                             reactions_to_disable=None):
    """Sets all of the input fluxes to the model to zero and then will
    initialize the input fluxes to the values specified in the_medium if
    it is a dict or will see if the model has a composition dict and use
    that to do the initialization.

    cobra_model: A cobra.Model object.


    the_medium: A string, or a dictionary.
    If a string then the initialize_growth_medium function expects that
    the_model has an attribute dictionary called media_compositions, which is a
    dictionary of dictionaries for various medium compositions.  Where a medium
    composition is a dictionary of external boundary reaction ids for the
    medium components and the external boundary fluxes for each medium
    component.


    external_boundary_compartment: None or a string.
    If not None then it specifies the compartment in which to disable all of
    the external systems boundaries.

    external_boundary_reactions: None or a list of external_boundaries that are
    to have their bounds reset.  This acts in conjunction with
    external_boundary_compartment.


    reaction_lower_bound: Float.  The default value to use for the lower
    bound for the boundary reactions.

    reaction_upper_bound: Float.  The default value to use for the upper
    bound for the boundary.

    irreversible: Boolean.  If the model is irreversible then the medium
    composition is taken as the upper bound

    reactions_to_disable: List of reactions for which the upper and lower
    bounds are disabled.  This is superceded by the contents of
    media_composition

    """
    # Zero all of the inputs to the model
    if hasattr(the_medium, 'keys'):
        medium_composition = the_medium
    else:
        if hasattr(cobra_model, 'media_compositions'):
            if the_medium in cobra_model.media_compositions:
                medium_composition = cobra_model.media_compositions[the_medium]
            else:
                raise Exception("%s is not in the model's media list" %
                                the_medium)
        else:
            raise Exception("the model doesn't have attribute "
                            "media_compositions and the medium is not a dict")
    if external_boundary_reactions is not None:
        if isinstance(external_boundary_reactions[0], str):
            external_boundary_reactions = map(cobra_model.reactions.get_by_id,
                                              external_boundary_reactions)
    elif external_boundary_compartment is None:
            warn("We are initializing the medium without first adjusting all"
                 "external boundary reactions")

    # Select the system_boundary reactions to reset
    if external_boundary_compartment is not None:
        _system_boundaries = dict([(x, x.get_compartments())
                                   for x in cobra_model.reactions
                                   if x.boundary == 'system_boundary'])
        [_system_boundaries.pop(k) for k, v in list(_system_boundaries.items())
         if len(v) == 1 and external_boundary_compartment not in v]
        if external_boundary_reactions is None:
            external_boundary_reactions = _system_boundaries.keys()
        else:
            external_boundary_reactions += _system_boundaries.keys()

    for the_reaction in external_boundary_reactions:
        the_reaction.lower_bound = reaction_lower_bound
        if the_reaction.upper_bound == 0:
            the_reaction.upper_bound = reaction_upper_bound
    # Disable specified reactions
    if reactions_to_disable is not None:
        if isinstance(reactions_to_disable[0], str):
            reactions_to_disable = map(cobra_model.reactions.get_by_id,
                                       reactions_to_disable)
        for the_reaction in reactions_to_disable:
            the_reaction.lower_bound = the_reaction.upper_bound = 0.

    # Update the model inputs based on the_medium
    for the_component in medium_composition.keys():
        the_reaction = cobra_model.reactions.get_by_id(the_component)
        if irreversible:
            the_reaction.upper_bound = medium_composition[the_component]
        else:
            the_reaction.lower_bound = medium_composition[the_component]


def convert_to_irreversible(cobra_model):
    """Split reversible reactions into two irreversible reactions

    These two reactions will proceed in opposite directions. This
    guarentees that all reactions in the model will only allow
    positive flux values, which is useful for some modeling problems.

    cobra_model: A Model object which will be modified in place.

    """
    reactions_to_add = []
    for reaction in cobra_model.reactions:
        # If a reaction is reverse only, the forward reaction (which
        # will be constrained to 0) will be left in the model.
        if reaction.lower_bound < 0:
            reverse_reaction = Reaction(reaction.id + "_reverse")
            reverse_reaction.lower_bound = 0
            reverse_reaction.upper_bound = reaction.lower_bound * -1
            reverse_reaction.objective_coefficient = \
                reaction.objective_coefficient * -1
            reaction.lower_bound = 0
            # Make the directions aware of each other
            reaction.reflection = reverse_reaction
            reverse_reaction.reflection = reaction
            reaction_dict = dict([(k, v*-1)
                                  for k, v in reaction._metabolites.items()])
            reverse_reaction.add_metabolites(reaction_dict)
            reverse_reaction._model = reaction._model
            reverse_reaction._genes = reaction._genes
            for gene in reaction._genes:
                gene._reaction.add(reverse_reaction)
            reverse_reaction._gene_reaction_rule = reaction._gene_reaction_rule
            reactions_to_add.append(reverse_reaction)
    cobra_model.add_reactions(reactions_to_add)


def revert_to_reversible(cobra_model, update_solution=True):
    """This function will convert a reversible model made by
    convert_to_irreversible into a reversible model.

    cobra_model: A cobra.Model which will be modified in place.

    """
    reverse_reactions = [x for x in cobra_model.reactions
                         if x.reflection is not None and
                         x.id.endswith('_reverse')]

    # If there are no reverse reactions, then there is nothing to do
    if len(reverse_reactions) == 0:
        return

    for reverse in reverse_reactions:
        forward = reverse.reflection
        forward.lower_bound = -reverse.upper_bound
        forward.reflection = None
    # Since the metabolites and genes are all still in
    # use we can do this faster removal step.  We can
    # probably speed things up here.
    cobra_model.remove_reactions(reverse_reactions)
    # fix the solution
    if update_solution and cobra_model.solution is not None and \
            cobra_model.solution.status != "NA":
        x_dict = cobra_model.solution.x_dict
        # Check if the model was optimized while it was still reversible. If so
        # then the solution does not need to be fixed.
        if reverse_reactions[0].id not in x_dict:
            return
        # Update x and x_dict to correct for reverse flux.
        for reverse in reverse_reactions:
            forward = reverse.reflection
            x_dict[forward.id] -= x_dict.pop(reverse.id)
        cobra_model.solution.x = [x_dict[r.id] for r in cobra_model.reactions]
