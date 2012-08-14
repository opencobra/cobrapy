#cobra.flux_analysis.objective.py
#functions for analyzing / creating objective functions
from ..core.Reaction import Reaction
from ..manipulation import initialize_growth_medium
import sys
if hasattr(sys, 'maxsize') and sys.maxsize > 2**32:
    try:
        from numpy import int64, int32
    except:
        int32 = int64 = int
else:
    int32 = int
    int64 = int
def assess_objective(cobra_model, the_objective=None,
                      objective_cutoff=0.001, growth_medium=None):
    """Assesses the ability of the model to produce all reactants in the_objective on
    an individual basis.  Returns True if the_objective can be realized to exceed
    objective_cutoff.  Otherwise, determines which components of the_objective are
    lagging and returns a dict of the components and their required and realized values.

    """
    cobra_model = cobra_model.copy()
    if growth_medium:
        initialize_growth_medium(cobra_model, growth_medium)

    #If the model cannot achieve the objective then check each component for failure
    #to be produced.
    if the_objective is None:
        #Assume a single objective reaction
        objective_reaction = [x for x in cobra_model.reactions if x.objective_coefficient != 0][0]
    elif hasattr(the_objective, 'id'):
        objective_reaction = cobra_model.reactions.get_by_id(the_objective.id) #need to get because we've copied the model
    elif isinstance(the_objective, str):
        objective_reaction = cobra_model.reactions.get_by_id(the_objective)
    else:
        #assume that it's an index
        objective_reaction = cobra_model.reactions[the_objective]
    cobra_model.optimize(new_objective = objective_reaction)
    #First see if the model can realize the objective
    if cobra_model.solution.f >= objective_cutoff:
        return {}
    components = objective_reaction.get_reactants()
    simulation_results = {}
    #TODO:  Speed this section up.  Possibly by modifying Model.optimize() to
    #use and updated S and reuse the basis.
    for the_component in objective_reaction.get_reactants():
        #add in a sink reaction for each component
        sink_reaction = Reaction('test_sink_' + the_component)
        #then simulate ability
        #then check it can exceed objective cutoff * component stoichiometric
        #coefficient.
        tmp_coefficient = objective_reaction.get_coefficient(the_component) 
        sink_reaction.add_metabolites(the_component, tmp_coefficient)
        sink_reaction.upper_bound = 1000
        cobra_model.add_reaction(sink_reaction)
        cobra_model.optimize(new_objective = sink_reaction.id)
        if objective_cutoff > cobra_model.solution.f:
            simulation_results.update({the_component:{'required':abs(tmp_coefficient*objective_cutoff), 'produced':cobra_model.solution.f}})
    return simulation_results

def update_objective(cobra_model, the_objectives):
    """Revised to take advantage of the new Reaction classes.

    cobra_model:  A cobra.Model

    the_objectives: A list or a dictionary.  If a list then
    a list of reactions for which the coefficient in the
    linear objective is set as 1.  If a dictionary then the
    key is the reaction and the value is the linear coefficient
    for the respective reaction.

    """
    #set the objective coefficients for each reaction to 0
    [setattr(x, 'objective_coefficient', 0.)
     for x in cobra_model.reactions]
    if isinstance(the_objectives, dict):
        for the_reaction, the_coefficient in the_objectives.iteritems():
            if isinstance(the_reaction, int):
                the_reaction = cobra_model.reactions[the_reaction]
            else:
                if hasattr(the_reaction, 'id'):
                    the_reaction = the_reaction.id
                the_reaction = cobra_model.reactions.get_by_id(the_reaction)
            the_reaction.objective_coefficient = the_coefficient
    else:
        #Allow for objectives to be constructed from multiple reactions
        if not isinstance(the_objectives, list) and \
               not isinstance(the_objectives, tuple):
            the_objectives = [the_objectives]
        for the_reaction in the_objectives:
            if isinstance(the_reaction, int):
                the_reaction = cobra_model.reactions[the_reaction]
            else:
                if hasattr(the_reaction, 'id'):
                    the_reaction = the_reaction.id
                the_reaction = cobra_model.reactions.get_by_id(the_reaction)
            the_reaction.objective_coefficient = 1.

