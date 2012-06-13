#cobra.flux_analysis.objective.py
#functions for analyzing / creating objective functions
from ..core.Reaction import Reaction
from ..manipulation import initialize_growth_medium
import sys
if sys.maxsize > 2**32:
    from numpy import int64, int32
else:
    int32 = int
    int64 = int
def assess_objective( cobra_model, the_objective = None,
                      objective_cutoff = 0.001, growth_medium = None ):
    """Assesses the ability of the model to produce all reactants in the_objective on
    an individual basis.  Returns True if the_objective can be realized to exceed
    objective_cutoff.  Otherwise, determines which components of the_objective are
    lagging and returns a dict of the components and their required and realized values.

    """
    raise Exception("This needs to be cleaned up.")
    cobra_model = cobra_model.copy()
    if growth_medium:
        initialize_growth_medium(cobra_model, growth_medium)

    #If the model cannot achieve the objective then check each component for failure
    #to be produced.
    if the_objective is None:
        objective_reaction = Reaction()
        objective_reaction.populate_from_cobra_model(cobra_model, cobra_model._objective_coefficients.nonzero()[0][0])
    elif (isinstance(the_objective, int) or  isinstance(the_objective, int64)
           or isinstance(the_objective, int32)  or isinstance(the_objective, str)):
        #TODO: When Model is updated to contain  Objects just
        #copy from there
        #objective_reaction  = deepcopy( the_model.cobra.reaction[ the_objective ] )
        objective_reaction = Reaction()
        objective_reaction.populate_from_cobra_model(cobra_model, the_objective)
    elif hasattr(the_objective, 'id'):
        objective_reaction = the_objective
    cobra_model.optimize(new_objective = objective_reaction)
    #First see if the model can realize the objective
    if cobra_model.solution.f >= objective_cutoff:
        return {}
    components = objective_reaction.get_reactants()
    component_indices = map(cobra_model.metabolites.index, components)
    simulation_results = {}
    #TODO:  Speed this section up.  Possibly by modifying Model.optimize() to
    #use and updated S and reuse the basis.
    for the_component, the_index in zip(components, component_indices):
        #add in a demand reaction for each component
        demand_reaction = Reaction('test_demand_' + the_component)
        #then simulate ability
        #then check it can exceed objective cutoff * component stoichiometric
        #coefficient.
        tmp_coeff = objective_reaction.get_coefficient(the_component) 
        demand_reaction.add_metabolites(the_component, tmp_coeff)
        demand_reaction.upper_bound = 1000
        cobra_model.add_reaction(demand_reaction)
        cobra_model.optimize(new_objective = demand_reaction.id)
        if objective_cutoff > cobra_model.solution.f:
            simulation_results.update({the_component:{'required':abs(tmp_coeff*objective_cutoff), 'produced':cobra_model.solution.f}})
    return simulation_results
def update_objective(cobra_model, the_objectives):
    """Revised to take advantage of the new Reaction classes.

    """
    from numpy import array
    #set the objective coefficients for each reaction to 0
    [setattr(x, 'objective_coefficient', 0.)
     for x in cobra_model.reactions]
    #Allow for objectives to be constructed from multiple reactions
    if not isinstance(the_objectives, list) and \
           not isinstance(the_objectives, tuple):
        the_objectives = [the_objectives]
    for the_objective in the_objectives:
        if not hasattr(the_objective,'id'):
            if isinstance(the_objective, str):
                the_objective = cobra_model.reactions.get_by_id(the_objective)
            elif isinstance(the_objective, int):
                the_objective = cobra_model.reactions[the_objective]
        the_objective.objective_coefficient = 1.
    #NOTE: _objective_coefficients is deprecated
    if hasattr(cobra_model, '_objective_coefficients'):
        cobra_model._objective_coefficients = array([x.objective_coefficient
                                                     for x in cobra_model.reactions])
