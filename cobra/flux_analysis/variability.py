from warnings import warn

from ..external.six import iteritems, string_types
from ..core.Metabolite import Metabolite
from ..solvers import solver_dict, get_solver_name

def flux_variability_analysis(cobra_model, reaction_list=None,
                              fraction_of_optimum=1.0, solver=None,
                              objective_sense="maximize", **solver_args):
    """Runs flux variability analysis to find max/min flux values

    cobra_model : :class:`~cobra.core.Model`:

    reaction_list : list of :class:`~cobra.core.Reaction`: or their id's
        The id's for which FVA should be run. If this is None, the bounds
        will be comptued for all reactions in the model.

    fraction_of_optimum : fraction of optimum which must be maintained.
        The original objective reaction is constrained to be greater than
        maximal_value * fraction_of_optimum

    solver : string of solver name
        If None is given, the default solver will be used.

    """
    if reaction_list is None and "the_reactions" in solver_args:
        reaction_list = solver_args.pop("the_reactions")
        from warnings import warn
        warn("the_reactions is deprecated. Please use reaction_list=")
    if reaction_list is None:
        reaction_list = cobra_model.reactions
    else:
        reaction_list = [cobra_model.reactions.get_by_id(i) if isinstance(i, string_types) else i for i in reaction_list]
    solver = solver_dict[get_solver_name() if solver is None else solver]
    lp = solver.create_problem(cobra_model)
    solver.solve_problem(lp, objective_sense=objective_sense)
    solution = solver.format_solution(lp, cobra_model)
    # set all objective coefficients to 0
    for i, r in enumerate(cobra_model.reactions):
        if r.objective_coefficient != 0:
            f = solution.x_dict[r.id]
            new_bounds = (f * fraction_of_optimum, f)
            solver.change_variable_bounds(lp, i, min(new_bounds), max(new_bounds))
            solver.change_variable_objective(lp, i, 0.)
    # perform fva
    fva_results = {}
    for r in reaction_list:
        i = cobra_model.reactions.index(r)
        fva_results[r.id] = {}
        solver.change_variable_objective(lp, i, 1.)
        solver.solve_problem(lp, objective_sense="maximize", **solver_args)
        fva_results[r.id]["maximum"] = solver.get_objective_value(lp)
        solver.solve_problem(lp, objective_sense="minimize", **solver_args)
        fva_results[r.id]["minimum"] = solver.get_objective_value(lp)
        # revert the problem to how it was before
        solver.change_variable_objective(lp, i, 0.)
    return fva_results

def flux_variability_analysis_legacy(cobra_model, fraction_of_optimum=1.,
                              objective_sense='maximize', the_reactions=None,
                              allow_loops=True, solver=None,
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
    from math import floor,ceil
    if solver is None:
        solver = get_solver_name()
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
        from copy import deepcopy
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
            for the_sense, the_description in iteritems(the_sense_dict):
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
            for k, v in iteritems(original_objectives):
                k.objective_coefficient = v
            objective_metabolite.remove_from_model()
    return variability_dict


def find_blocked_reactions(cobra_model, the_reactions=None, allow_loops=True,
                            solver=None, the_problem='return',
                           tolerance_optimality=1e-9,
                           open_exchanges=False, **kwargs):
    """Finds reactions that cannot carry a flux with the current
    exchange reaction settings for cobra_model, using flux variability
    analysis.
    
    """
    if solver is None:
        solver = get_solver_name()
    warn('This needs to be updated to deal with external boundaries')
    cobra_model = cobra_model.copy()
    blocked_reactions = []
    if not the_reactions:
        the_reactions = cobra_model.reactions
    if open_exchanges:
        warn('DEPRECATED: Move to using the Reaction.boundary attribute')
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


