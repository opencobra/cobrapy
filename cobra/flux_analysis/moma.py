from __future__ import print_function
#cobra.flux_analysis.moma.py: Runs the minimization of metabolic
#adjustment method described in Segre et al 2002 PNAS 99(23): 15112-7
from os import name as __name
from sys import modules as __modules
from warnings import warn
if __name == 'java':
    raise Exception("%s is not yet supported on jython"%__modules[__name__])
from copy import deepcopy
from time import time
from math import ceil, floor

#The next four imports need to be dealt with to obtain jython compatibilty
from numpy import array, hstack, vstack, matrix, sum
from scipy.sparse import eye, lil_matrix, dok_matrix
from scipy.sparse import hstack as s_hstack
from scipy.sparse import vstack as s_vstack

from ..core import Reaction, Metabolite
from ..manipulation import initialize_growth_medium, delete_model_genes
from ..manipulation.modify import convert_to_irreversible
from ..external.six import iteritems

from warnings import warn

def moma(wt_model, mutant_model, objective_sense='maximize', solver=None,
         tolerance_optimality=1e-8, tolerance_feasibility=1e-8,
         minimize_norm=False, the_problem='return', lp_method=0,
         combined_model=None, norm_type='euclidean'):
    """Runs the minimization of metabolic adjustment method described in
    Segre et al 2002 PNAS 99(23): 15112-7.

    wt_model: A cobra.Model object

    mutant_model: A cobra.Model object with different reaction bounds vs wt_model.
    To simulate deletions

    objective_sense: 'maximize' or 'minimize'

    solver: 'gurobi', 'cplex', or 'glpk'.  Note: glpk cannot be used with
    norm_type 'euclidean'

    tolerance_optimality: Solver tolerance for optimality.

    tolerance_feasibility: Solver tolerance for feasibility.

    the_problem: None or a problem object for the specific solver that can be
    used to hot start the next solution.

    lp_method: The method to use for solving the problem.  Depends on the solver.  See
    the cobra.flux_analysis.solvers.py file for more info.
        For norm_type == 'euclidean':
            the primal simplex works best for the test model (gurobi: lp_method=0, cplex: lp_method=1)
    
    combined_model: an output from moma that represents the combined optimization
    to be solved.  when this is not none.  only assume that bounds have changed
    for the mutant or wild-type.  This saves 0.2 seconds in stacking matrices.


    NOTE: Current function makes too many assumptions about the structures of the models


    """
    if solver is None:
        if norm_type == "euclidiean":
            solver = get_solver_name(qp=True)
        else:
            solver = get_solver_name()  # linear is not even implemented yet
    if combined_model is not None or the_problem not in ['return']:
        warn("moma currently does not support reusing models or problems. " +\
             "continuing without them")
        combined_model = None
        the_problem = 'return'
    if solver.lower() == 'cplex' and lp_method == 0:
        #print 'for moma, solver method 0 is very slow for cplex. changing to method 1'
        lp_method = 1
    if solver.lower() == 'glpk' and norm_type == 'euclidean':
        try:
            from gurobipy import Model
            solver = 'gurobi'
            warn("GLPK can't solve quadratic problems like MOMA.  Switched solver to %s"%solver)
        except:
            warn("GLPK can't solve quadratic problems like MOMA.  Switching to linear MOMA")

    if norm_type == 'euclidean':
        #Reusing the basis can get the solver stuck.
        reuse_basis = False
    if combined_model and combined_model.norm_type != norm_type:
        print('Cannot use combined_model.norm_type = %s with user-specified norm type'%(combined_model.norm_type,
                                                                                        norm_type))
        print('Defaulting to user-specified norm_type')
        combined_model = None



    number_of_reactions_in_common = len(set([x.id for x in wt_model.reactions]).intersection([x.id for x in mutant_model.reactions]))
    number_of_reactions = len(wt_model.reactions) + len(mutant_model.reactions)

    #Get the optimal wt objective value and adjust based on optimality tolerances
    wt_model.optimize(solver=solver)
    wt_optimal = deepcopy(wt_model.solution.f)
    if objective_sense == 'maximize':
        wt_optimal = floor(wt_optimal/tolerance_optimality)*tolerance_optimality
    else:
        wt_optimal = ceil(wt_optimal/tolerance_optimality)*tolerance_optimality
    if norm_type == 'linear':
        raise Exception('linear MOMA is not currently implmented')
        quadratic_component = None
    if minimize_norm:
        raise Exception('minimize_norm is not currently implemented')
        #just worry about the flux distribution and not the objective from the wt
        combined_model = mutant_model.copy()
        #implement this: combined_model.reactions[:].objective_coefficients = -wt_solution.x_dict
    else:
        #Construct a problem that attempts to maximize the objective in the WT model while
        #solving the quadratic problem.  This new problem is constructed to try to find
        #a solution for the WT model that lies close to the mutant model.  There are
        #often multiple equivalent solutions with M matrices and the one returned
        #by a simple cobra_model.optimize call may be too far from the mutant.
        #This only needs to be adjusted if we update mutant_model._S after deleting reactions
        if not combined_model:
            #Collect the set of wt reactions contributing to the objective.
            objective_reaction_coefficient_dict = dict([(x.id, x.objective_coefficient)
                                                        for x in wt_model.reactions
                                                        if x.objective_coefficient])
            
            
            combined_model = construct_difference_model(wt_model, mutant_model, norm_type)
            #Add in the virtual objective metabolite to constrain the wt_model to the space where
            #the objective was maximal
            objective_metabolite = Metabolite('wt_optimal')
            objective_metabolite._bound = wt_optimal
            if objective_sense == 'maximize':
                objective_metabolite._constraint_sense = 'G'
            else:
                objective_metabolite._constraint_sense = 'L'

            #TODO: this couples the wt_model objective reaction to the virtual metabolite
            #Currently, assumes a single objective reaction; however, this may be extended
            [combined_model.reactions.get_by_id(k).add_metabolites({objective_metabolite: v})
             for k, v in objective_reaction_coefficient_dict.items()]
                
            
            



        if norm_type == 'euclidean':
            #Makes assumptions about the structure of combined model
            quadratic_component = s_vstack((lil_matrix((number_of_reactions, number_of_reactions + number_of_reactions_in_common )),
                                            s_hstack((lil_matrix((number_of_reactions_in_common, number_of_reactions)),
                                                      eye(number_of_reactions_in_common,number_of_reactions_in_common)))))
    
        elif norm_type == 'linear':
            quadratic_component = None

    combined_model.norm_type = norm_type
    cobra_model = combined_model

    the_problem = combined_model.optimize(objective_sense='minimize',
                                         quadratic_component=quadratic_component,
                                         solver=solver,
                                         tolerance_optimality=tolerance_optimality,
                                         tolerance_feasibility=tolerance_feasibility,
                                         lp_method=lp_method, reuse_basis=reuse_basis)

    if combined_model.solution.status != 'optimal':
        warn('optimal moma solution not found: solver status %s'%combined_model.solution.status +\
             ' returning the problem, the_combined model, and the quadratic component for trouble shooting')
        return(the_problem, combined_model, quadratic_component)
             
    solution = combined_model.solution
    mutant_dict = {}
    #Might be faster to quey based on mutant_model.reactions with the 'mutant_' prefix added
    _reaction_list = [x for x in combined_model.reactions if x.id.startswith('mutant_')]
    mutant_f = sum([mutant_model.reactions.get_by_id(x.id[len('mutant_'):]).objective_coefficient *
                    x.x for x in _reaction_list])
    mutant_dict['objective_value'] = mutant_f
    wild_type_flux_total = sum([abs(solution.x_dict[x.id]) for x in wt_model.reactions])
    mutant_flux_total = sum(abs(x.x) for x in _reaction_list)
    #Need to use the new solution as there are multiple ways to achieve an optimal solution in
    #simulations with M matrices.
    mutant_dict['status'] = solution.status
    #TODO: Deal with maximize / minimize issues for a reversible model that's been converted to irreversible
    mutant_dict['flux_difference'] = flux_difference = sum([(solution.x_dict[x.id[len('mutant_'):]]
                                                             - x.x)**2 for x in _reaction_list])
    mutant_dict['the_problem'] = the_problem
    mutant_dict['combined_model'] = combined_model
    
    del wt_model, mutant_model, quadratic_component, solution
    return(mutant_dict)


def construct_difference_model(model_1, model_2, norm_type='euclidean'):
    """Combine two models into a larger model that is designed to calculate differences
    between the models

    """
    #Get index mappings
    common_dict = {}
    #Using copies of the models so things are modified above
    combined_model = model_1 = model_1.copy()
    model_2 = model_2.copy()
    for reaction_1 in model_1.reactions:
        try:
            reaction_2 = model_2.reactions.get_by_id(reaction_1.id)
            common_dict[reaction_1] = reaction_2
        except:
            continue
            
    #Add a prefix in front of the mutant_model metabolites and reactions to prevent
    #name collisions in DictList
    for the_dict_list in [model_2.metabolites,
                          model_2.reactions]:
        [setattr(x, 'id', 'mutant_%s'%x.id)
         for x in the_dict_list]
        the_dict_list._generate_index() #Update the DictList.dicts

    
    combined_model.add_reactions(model_2.reactions)
    [setattr(x, 'objective_coefficient', 0.)
     for x in combined_model.reactions]
    #Add in the difference reactions.  The mutant reactions and metabolites are already added.
    #This must be a list to maintain the correct order when adding the difference_metabolites
    difference_reactions = [] #Add the difference reactions at the end to speed things up
    difference_metabolites = []
    for reaction_1, reaction_2 in iteritems(common_dict):
        reaction_1._difference_partner = reaction_2
        reaction_2._difference_partner = reaction_1
        difference_reaction = Reaction('difference_%s'%reaction_1.id)
        difference_reactions.append(difference_reaction)
        difference_reaction.upper_bound = 100000
        difference_reaction.lower_bound = -1* difference_reaction.upper_bound
        difference_metabolite = Metabolite('difference_%s'%reaction_1.id)
        difference_metabolites.append(difference_metabolite)
        if norm_type == 'linear':
            difference_metabolite._constraint_sense = 'G'
        reaction_1.add_metabolites({difference_metabolite: -1.}, add_to_container_model=False)
        reaction_2.add_metabolites({difference_metabolite: 1.}, add_to_container_model=False)
        difference_reaction.add_metabolites({difference_metabolite: 1.}, add_to_container_model=False)

    combined_model.add_metabolites(difference_metabolites)
    combined_model.add_reactions(difference_reactions)
    return(combined_model)

