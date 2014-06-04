from time import time
from warnings import warn
from copy import deepcopy

from ..external.six import string_types, iteritems

from ..manipulation import delete_model_genes, undelete_model_genes
from ..manipulation.delete import find_gene_knockout_reactions
from ..solvers import solver_dict, get_solver_name


nan = float('nan')


try:
    from .moma import moma
except Exception as e:
    def moma(**kwargs):
        warn("moma is currently not functional")


def single_deletion(cobra_model, element_list=None,
                    method='fba', element_type='gene', solver=None):
    """Wrapper for single_gene_deletion and the single_reaction_deletion
    functions

    
    cobra_model: a cobra.Model object

    element_list: Is None or a list of elements (genes or reactions) to
    delete.

    method: 'fba' or 'moma'

    the_problem: Is None, 'return', or an LP model object for the solver.

    element_type: 'gene' or 'reaction'

    solver: 'glpk', 'gurobi', or 'cplex'.

    error_reporting: None or True to disable or enable printing errors encountered
    when trying to find the optimal solution.

    discard_problems: Boolean.  If True do not save problems.  This will
    help with memory issues related to gurobi.
    .. warning:: This is deprecated.

    Returns a list of dictionaries: growth_rate_dict, solution_status_dict,
    problem_dict where the key corresponds to each element in element_list.

    """
    if solver is None:
        solver = get_solver_name() if method == "fba" else get_solver_name(qp=True)
    # fast versions of functions
    if method == "fba":
        if element_type == "gene":
            return single_gene_deletion_fba(cobra_model, element_list,
                                             solver=solver)
        elif element_type == "reaction":
            return single_reaction_deletion_fba(cobra_model, element_list,
                                                 solver=solver)
    if element_type == 'gene':
        the_solution = single_gene_deletion(cobra_model, element_list,
                                    method=method, solver=solver)

    else:
        the_solution = single_reaction_deletion(cobra_model, element_list,
                                        method=method, solver=solver)
    return the_solution


def single_reaction_deletion_fba(cobra_model, reaction_list=None, solver=None):
    solver = solver_dict[get_solver_name() if solver is None else solver]
    lp = solver.create_problem(cobra_model)
    growth_rate_dict = {}
    status_dict = {}
    if reaction_list is None:
        reaction_list = cobra_model.reactions
    else:
        reaction_list = [cobra_model.reactions.get_by_id(i) \
                         if isinstance(i, string_types) else i \
                         for i in reaction_list]
    for reaction in reaction_list:
        old_bounds = (reaction.lower_bound, reaction.upper_bound)
        index = cobra_model.reactions.index(reaction)
        solver.change_variable_bounds(lp, index, 0., 0.)
        solver.solve_problem(lp)
        status = solver.get_status(lp)
        status_dict[reaction.id] = status
        growth_rate_dict[reaction.id] = solver.get_objective_value(lp) if status == "optimal" else 0.
        # reset the problem
        solver.change_variable_bounds(lp, index, old_bounds[0], old_bounds[1])
    return(growth_rate_dict, status_dict)

def single_gene_deletion_fba(cobra_model, gene_list=None, solver=None):
    solver = solver_dict[get_solver_name() if solver is None else solver]
    lp = solver.create_problem(cobra_model)
    growth_rate_dict = {}
    status_dict = {}
    if gene_list is None:
        gene_list = cobra_model.genes
    else:
        gene_list = [cobra_model.genes.get_by_id(i) \
                     if isinstance(i, string_types) else i for i in gene_list]
    for gene in gene_list:
        old_bounds = {}
        for reaction in find_gene_knockout_reactions(cobra_model, [gene]):
            index = cobra_model.reactions.index(reaction)
            old_bounds[index] = (reaction.lower_bound, reaction.upper_bound)
            solver.change_variable_bounds(lp, index, 0., 0.)
        solver.solve_problem(lp)
        status = solver.get_status(lp)
        status_dict[gene.id] = status
        growth_rate_dict[gene.id] = solver.get_objective_value(lp)# if status == "optimal" else 0.
        # reset the problem
        for index, bounds in iteritems(old_bounds):
            solver.change_variable_bounds(lp, index, bounds[0], bounds[1])
    return(growth_rate_dict, status_dict)

def single_reaction_deletion(cobra_model, element_list=None,
                             method='fba', the_problem='return',
                             solver=None, error_reporting=None,
                             discard_problems=True):
    """Performs optimization simulations to realize the objective defined
    from cobra_model.reactions[:].objective_coefficients after deleting each reaction
    from the model.
    
    cobra_model: a cobra.Model object

    element_list: Is None or a list of cobra.Reactions in cobra_model to disable.
    If None then disable each reaction in cobra_model.reactions and optimize for the
    objective function defined from cobra_model.reactions[:].objective_coefficients.

    method: 'fba' is the only option at the moment.

    the_problem: Is None, 'reuse', or an LP model object for the solver.

    solver: 'glpk', 'gurobi', or 'cplex'.

    discard_problems: Boolean.  If True do not save problems.  This will
    help with memory issues related to gurobi.
    
    Returns a list of dictionaries: growth_rate_dict, solution_status_dict,
    problem_dict where the key corresponds to each reaction in reaction_list.

    """
    if solver is None:
        solver = get_solver_name() if method == "fba" else get_solver_name(qp=True)
    #element_list so we can merge single_reaction_deletion and single_gene_deletion

    wt_model = cobra_model.copy() #Original wild-type (wt) model.
    wt_model.id = 'Wild-Type'
    #MOMA constructs combined quadratic models thus we cannot reuse a model
    #generated by the cobra_model.optimize call
    if method.lower() == 'moma':
        the_problem = 'return'
        mutant_model = wt_model.copy() #Need a second model for moma
    else:
        mutant_model = cobra_model
    discard_problems = False
    if the_problem:
        the_problem = 'return'
        discard_problems = True
    the_problem = wt_model.optimize(the_problem=the_problem, solver=solver,
                                       error_reporting=error_reporting)
    wt_f = wt_model.solution.f
    wt_status = wt_model.solution.status
    wt_x = deepcopy(wt_model.solution.x)
    wt_x_dict = deepcopy(wt_model.solution.x_dict)

    wt_problem = the_problem
    if element_list is None:
        element_list = mutant_model.reactions
    elif not hasattr(element_list[0], 'id'):
        element_list = map(mutant_model.reactions.get_by_id, element_list)
    else:
        if mutant_model is not cobra_model:
            element_list = [x.id for x in element_list]
            element_list = map(mutant_model.genes.get_by_id, element_list)

    growth_rate_dict = {}
    solution_status_dict = {}
    problem_dict = {}
    combined_model = None
    for the_element in element_list:
        #delete the gene
        #if the deletion alters the bounds then run simulation
        old_lower_bound, old_upper_bound = map(float, [the_element.lower_bound,
                                                       the_element.upper_bound])
        mutant_model.id = the_element.id
        if old_lower_bound != 0 or old_upper_bound != 0:
            the_element.lower_bound = 0
            the_element.upper_bound = 0
            if method.lower() == 'fba':
                the_problem = mutant_model.optimize(the_problem=wt_problem,
                                                    solver=solver,
                                                    error_reporting=error_reporting)
                growth_rate_dict[the_element] = mutant_model.solution.f
                solution_status_dict[the_element] = mutant_model.solution.status
            elif method.lower() == 'moma':
                try:
                    #TODO: Need to figure out why reusing the problem and the combined_model do not
                    #speed things up here.
                    moma_solution = moma(wt_model, mutant_model, solver=solver, the_problem=the_problem,
                                         combined_model=combined_model)
                    the_problem = moma_solution.pop('the_problem')
                    growth_rate_dict[the_element] = float(moma_solution.pop('objective_value'))
                    solution_status_dict[the_element] = moma_solution.pop('status')
                    combined_model = moma_solution.pop('combined_model')
                except:
                    growth_rate_dict[the_element] = nan
                    the_problem = None
                    solution_status_dict[the_element] = 'failed'

            if discard_problems:
                problem_dict[the_element] = 'discarded'
            else:
                problem_dict[the_element] = the_problem
            if not the_problem:
                the_problem = wt_problem
            #reset the model
            the_element.lower_bound = old_lower_bound
            the_element.upper_bound = old_upper_bound
        #else just use the wt_f and x
        else:
            if discard_problems:
                problem_dict[the_element] = 'discarded'
            else:
                problem_dict[the_element] = wt_problem
            growth_rate_dict[the_element] = wt_f
            solution_status_dict[the_element] = wt_status
    
    
    return(growth_rate_dict, solution_status_dict, problem_dict)

def single_gene_deletion(cobra_model, element_list=None,
                         method='fba', the_problem='reuse', solver=None,
                         error_reporting=None):
    """Performs optimization simulations to realize the objective defined
    from cobra_model.reactions[:].objective_coefficients after deleting each gene in
    gene_list from the model.
    
    cobra_model: a cobra.Model object

    element_list: Is None or a list of genes to delete.  If None then
    disable each reaction associated with each gene in cobra_model.genes.

    method: 'fba' or 'moma'

    the_problem: Is None or 'reuse'.

    solver: 'glpk', 'gurobi', or 'cplex'.

    Returns a list of dictionaries: growth_rate_dict, solution_status_dict,
    problem_dict where the key corresponds to each reaction in reaction_list.

    TODO: Add in a section that allows copying and collection of problem for
    debugging purposes.

    """
    if solver is None:
        solver = get_solver_name() if method == "fba" else get_solver_name(qp=True)
    wt_model = cobra_model.copy() #Original wild-type (wt) model.
    wt_model.id = 'Wild-Type'
    #MOMA constructs combined quadratic models thus we cannot reuse a model
    #generated by the cobra_model.optimize call
    if method.lower() == 'moma':
        the_problem = 'return'
        mutant_model = wt_model.copy() #Need a second model for moma
    else:
        mutant_model = cobra_model
    discard_problems = False
    if the_problem:
        the_problem = 'return'
        discard_problems = True

    solver_object = solver_dict[solver]
    the_problem = solver_object.create_problem(wt_model)
    solver_object.solve_problem(the_problem)
    solution = solver_object.format_solution(the_problem, wt_model)
    wt_f = solution.f
    wt_status = solution.status
    wt_x = deepcopy(solution.x)
    wt_x_dict = deepcopy(solution.x_dict)

    if element_list is None:
        element_list = mutant_model.genes
    elif not hasattr(element_list[0], 'id'):
        element_list = map(mutant_model.genes.get_by_id, element_list)
    else:
        if mutant_model is not cobra_model:
            element_list = [x.id for x in element_list]
            element_list = map(mutant_model.genes.get_by_id, element_list)
    wt_problem = the_problem

    growth_rate_dict = {}
    solution_status_dict = {}
    problem_dict = {}
    combined_model = None

    for the_element in element_list:
        #delete the gene
        #if the deletion alters the bounds then run simulation
        delete_model_genes(mutant_model, the_element)
        mutant_model.id = the_element.id
        if mutant_model._trimmed:
            if method.lower() == 'fba':
                the_problem = mutant_model.optimize(the_problem=wt_problem,
                                                    solver=solver,
                                                    error_reporting=error_reporting)
                growth_rate_dict[the_element.id] = mutant_model.solution.f
                solution_status_dict[the_element.id] = mutant_model.solution.status
            elif method.lower() == 'moma':
                try:
                    #TODO: Need to figure out why reusing the problem and the combined_model do not
                    #speed things up here.
                    moma_solution = moma(wt_model, mutant_model, solver=solver, the_problem=the_problem,
                                         combined_model=combined_model)
                    the_problem = moma_solution.pop('the_problem')
                    growth_rate_dict[the_element.id] = float(moma_solution.pop('objective_value'))
                    solution_status_dict[the_element.id] = moma_solution.pop('status')
                    combined_model = moma_solution.pop('combined_model')
                except:
                    growth_rate_dict[the_element.id] = nan
                    the_problem = None
                    solution_status_dict[the_element.id] = 'failed'

            if discard_problems:
                problem_dict[the_element.id] = 'discarded'
            else:
                problem_dict[the_element.id] = the_problem
            if not the_problem:
                the_problem = wt_problem
            #reset the model
            undelete_model_genes(mutant_model)
        #else just use the wt_f and x
        else:
            if discard_problems:
                problem_dict[the_element.id] = 'discarded'
            else:
                problem_dict[the_element.id] = wt_problem
            growth_rate_dict[the_element.id] = wt_f
            solution_status_dict[the_element.id] = wt_status
    
    
    return(growth_rate_dict, solution_status_dict, problem_dict)

