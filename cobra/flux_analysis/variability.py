from __future__ import with_statement
#cobra.flux_analysis.variablity.py
#runs flux variablity analysis on a Model object.
from math import floor,ceil
from copy import deepcopy
from ..core.Metabolite import Metabolite
try:
    from ..external.ppmap import ppmap
    __parallel_mode_available = True
except:
    __parallel_mode_available = False
#TODO: Add in a ppmap section for running in parallel
def flux_variability_analysis_wrapper(keywords):
    """Provides an interface to call flux_variability_analysis from ppmap
    
    """
    try:
        from cPickle import dump
    except:
        from pickle import dump
    import sys
    #Need to do a top level import because of how the parallel processes are called
    sys.path.insert(0, "../..")
    from cobra.flux_analysis.variability import flux_variability_analysis
    sys.path.pop(0)
    results_dict = {}
    new_objective = keywords.pop('new_objective')
    output_directory = None
    if 'output_directory' in keywords:
        output_directory = keywords.pop('output_directory')
    if not hasattr(new_objective, '__iter__'):
        new_objective = [new_objective]
    for the_objective in new_objective:
        the_result = results_dict[the_objective] = flux_variability_analysis(**keywords)
        if output_directory:
            with open('%s%s.pickle'%(output_directory,
                                     the_objective), 'w') as out_file:
                dump(the_result, out_file)
    if len(new_objective) == 1:
        return the_result
    else:
        return results_dict

def flux_variability_analysis(cobra_model, fraction_of_optimum=1.,
                              objective_sense='maximize', the_reactions=None,
                              allow_loops=True, solver='glpk',
                              the_problem='return', tolerance_optimality=1e-6,
                              tolerance_feasibility=1e-6, tolerance_barrier=1e-8,
                              lp_method=1, lp_parallel=0, new_objective=None,
                              relax_b=None, error_reporting=None,
                              number_of_processes=1, copy_model=True):
    """Runs flux variability analysis on a cobra.Model object

    cobra_model: a Model object

    fraction_of_optimum: fraction of the optimal solution that must be realized

    the_reactions: list of reactions to run FVA on.  if None then run on all
    reactions in the Model

    allow_loops:  Not Implemented.  If false then run the simulations with the
    loop law method to remove loops.

    the_problem: If 'return' or an LP model object for the specified solver then
    the optimizations will be sped up by attempting to use a previous solution
    as a starting point to optimize the current problem.  Can reduce
    simulation time by over an order of magnitude.

    solver: 'glpk', 'gurobi', or 'cplex'.

    the_problem: a problem object for the corresponding solver, 'return', or
    a float representing the wt_solution

    number_of_processes: If greater than 1 then this function will attempt
    to parallelize the problem.  NOTE: Currently not functional


    returns a dictionary: {reaction.id: {'maximum': float, 'minimum': float}}

    TODO: update how Metabolite._bound is handled so we can set a range instead
    of just a single value.  This will be done in cobra.flux_analysis.solvers.
    
    """
    #Need to copy the model because we're updating reactions.  However,
    #we can always just remove them.
    if isinstance(the_problem, float):
        wt_solution = the_problem
        the_problem='return'
    else:
        wt_model = cobra_model
        wt_model.optimize(solver=solver,objective_sense='maximize',
                          tolerance_optimality=tolerance_optimality,
                          tolerance_feasibility=tolerance_feasibility,
                          tolerance_barrier=tolerance_barrier,
                          lp_method=lp_method, lp_parallel=lp_parallel,
                          the_problem=the_problem, new_objective=new_objective)
        wt_solution = wt_model.solution.f
    if copy_model:
        cobra_model = cobra_model.copy()
    if not the_reactions:
        the_reactions = cobra_model.reactions
    else:
        if hasattr(the_reactions[0], 'id'):
            #Because cobra_model = cobra_model.copy() any cobra.Reactions
            #from the input won't point to cobra_model
            the_reactions = [x.id for x in the_reactions]
        #
        the_reactions = map(cobra_model.reactions.get_by_id, the_reactions)
    #If parallel mode is called for then give it a try
    if number_of_processes > 1 and __parallel_mode_available:
        the_problem = wt_solution #Solver objects are not thread safe
        the_reactions = [x.id for x in the_reactions]
        parameter_dict = dict([(x, eval(x))
                               for x in ['fraction_of_optimum',
                                         'objective_sense',
                                         'allow_loops',
                                         'solver',
                                         'the_problem',
                                         'tolerance_optimality',
                                         'tolerance_feasibility',
                                         'tolerance_barrier',
                                         'lp_method',
                                         'lp_parallel',
                                         'new_objective',
                                         'relax_b',
                                         'error_reporting']])
        #Might need to deepcopy when threading
        parameter_dict['cobra_model'] = cobra_model.copy()
        parameter_dict['copy_model'] = False
        parameter_list = []
        number_of_reactions = len(the_reactions)
        if number_of_reactions < number_of_processes:
            number_of_processes = number_of_reactions
        reactions_per_process = int(floor(number_of_reactions/number_of_processes))
        for i in range(number_of_processes):
            tmp_parameters = deepcopy(parameter_dict)
            tmp_parameters['the_reactions'] = the_reactions[:reactions_per_process]
            parameter_list.append(tmp_parameters)
            the_reactions = the_reactions[reactions_per_process:]
        #The remainder goes with the last processes
        parameter_list[-1]['the_reactions'] += the_reactions
        variability_dict = {}
        pp_pointer = ppmap(number_of_processes,
                           flux_variability_analysis_wrapper, parameter_list)

        [variability_dict.update(x) for x in list(pp_pointer)]
    else:
        #Basically, add a virtual metabolite that reflects the
        #objective coefficeints and the solution
        objective_metabolite = Metabolite('objective')
        if not copy_model:
            original_objectives = dict([(k, float(k.objective_coefficient))
                                        for k in cobra_model.reactions])
                                        
        [x.add_metabolites({objective_metabolite: x.objective_coefficient})
         for x in cobra_model.reactions if x.objective_coefficient != 0]

        #TODO: Kick back an error if cobra_model.solution.status is not optimal
        if objective_sense == 'maximize':
            objective_cutoff = floor(wt_solution/tolerance_optimality)*\
                               tolerance_optimality*fraction_of_optimum
            objective_metabolite._constraint_sense = 'G'
        else:
            objective_cutoff = ceil(wt_solution/tolerance_optimality)*\
                               tolerance_optimality*fraction_of_optimum
            objective_metabolite._constraint_sense = 'L'
        objective_metabolite._bound = objective_cutoff
        #If objective_metabolite._model is None then we should cycle through
        #each reaction as the initial objective.
        if the_problem:
            the_problem = cobra_model.optimize(solver=solver,
                                               objective_sense='maximize',
                                               tolerance_optimality=tolerance_optimality,
                                               tolerance_feasibility=tolerance_feasibility,
                                               tolerance_barrier=tolerance_barrier,
                                               lp_method=lp_method, lp_parallel=lp_parallel,
                                               the_problem='return', relax_b=relax_b,
                                               error_reporting=error_reporting)
        variability_dict = {}
        the_sense_dict = {'maximize': 'maximum',
                          'minimize': 'minimum'}
        basic_problem = the_problem
        #Adding in solver-specific code could improve the speed substantially.
        for the_reaction in the_reactions:
            tmp_dict = {}
            the_problem = basic_problem
            for the_sense, the_description in the_sense_dict.iteritems():
                the_problem = cobra_model.optimize(solver=solver,
                                                   new_objective=the_reaction,
                                                   objective_sense=the_sense,
                                                   tolerance_optimality=tolerance_optimality,
                                                   tolerance_feasibility=tolerance_feasibility,
                                                   tolerance_barrier=tolerance_barrier,
                                                   lp_method=lp_method,
                                                   lp_parallel=lp_parallel,
                                                   the_problem=the_problem,
                                                   error_reporting=error_reporting,
                                                   update_problem_reaction_bounds=False)
                tmp_dict[the_description] = cobra_model.solution.f
            variability_dict[the_reaction.id] = tmp_dict
        if not copy_model:
            [setattr(k, 'objective_coefficient', v)
             for k, v in original_objectives.iteritems()]
            objective_metabolite.remove_from_model()
    return variability_dict


def find_blocked_reactions(cobra_model, the_reactions=None, allow_loops=True,
                            solver='glpk', the_problem='return',
                           tolerance_optimality=1e-9,
                           open_exchanges=False, **kwargs):
    """Finds reactions that cannot carry a flux with the current
    exchange reaction settings for cobra_model, using flux variability
    analysis.
    
    """
    print 'This needs to be updated to deal with external boundaries'
    cobra_model = cobra_model.copy()
    blocked_reactions = []
    if not the_reactions:
        the_reactions = cobra_model.reactions
    if open_exchanges:
        print 'DEPRECATED: Move to using the Reaction.boundary attribute'
        exchange_reactions = [x for x in cobra_model.reactions
                              if x.startswith('EX')]
        for the_reaction in exchange_reactions:
            if the_reaction.lower_bound >= 0:
                the_reaction.lower_bound = -1000
            if the_reaction.upper_bound >= 0:
                the_reaction.upper_bound = 1000
    flux_span_dict = flux_variability_analysis(cobra_model,
                                               fraction_of_optimum = 0.,
                                               the_reactions = the_reactions,
                                               allow_loops = allow_loops,
                                               solver = solver,
                                               the_problem = the_problem,
                                               tolerance_optimality = tolerance_optimality,
                                               **kwargs)
    blocked_reactions = [k for k, v in flux_span_dict.items()\
                          if max(map(abs,v.values())) < tolerance_optimality]
    return(blocked_reactions)


def run_fva(solver='cplex'):
    """
    For debugging purposes only
    
    """
    from test import salmonella_pickle
    try:
        from cPickle import load
    except:
        from pickle import load
    with open(salmonella_pickle) as in_file:
        cobra_model = load(in_file)
    fva_out =  flux_variability_analysis(cobra_model,
                                         the_reactions=cobra_model.reactions,#[100:140],
                                         solver=solver,number_of_processes=1)

                             

