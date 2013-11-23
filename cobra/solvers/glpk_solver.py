##cobra.solvers.glpk_solver
#This script provides wrappers for libglpk-java 1.0.22 and pyglpk 0.3
from warnings import warn
from copy import deepcopy
###solver specific parameters
from .parameters import parameter_mappings, parameter_defaults, \
     default_objective_sense

from ..core.Solution import Solution
from ..flux_analysis.objective import update_objective
solver_name = 'glpk'


from glpk import LPX as GLPK
__solver_class = GLPK

variable_kind_dict = {'continuous': float, 'integer': int}
status_dict = {'opt': 'optimal', 'nofeas': 'infeasible', 'unbnd': 'unbounded'}
sense_dict = {'E': 'E', 'L': 'L', 'G': 'G'}
parameter_mappings = parameter_mappings[solver_name]
parameter_defaults = parameter_defaults[solver_name]

def get_status(lp):
    status = lp.status
    if status in status_dict:
        status = status_dict[status]
    else:
        status = 'failed'
    return status

def get_objective_value(lp):
    return lp.obj.value

def format_solution(lp, cobra_model, **kwargs):
    status = get_status(lp)
    if status == 'optimal':
        sol = Solution(lb.obj.value, status=status)
        sol.x = [float(c.primal) for c in lp.cols]
        sol.x_dict = {c.name: c.primal for c in lp.cols}

        # return the duals as well as the primals for LPs
        if lp.kind == float:
            sol.y = [float(c.dual) for c in lp.rows]
            y_dict = {c.name: c.dual for c in lp.rows}
        return sol

    return Solution(None, status=status)

def set_parameter(lp, parameter_name, parameter_value):
    """with pyglpk the parameters are set during the solve phase, with
    the exception of objective sense.
    
    """
    if parameter_name == 'objective_sense':
        if parameter_value.lower() == 'maximize':
            lp.obj.maximize = True
        elif parameter_value.lower() == 'minimize':
            lp.obj.maximize = False
        else:
            raise ValueError("objective_sense should be 'maximize' or 'minimize'")
    else:
        warn("py glpk solver parameters are set during solve_problem")


def create_problem(cobra_model,  **kwargs):
    """Solver-specific method for constructing a solver problem from
    a cobra.Model.  This can be tuned for performance using kwargs


    """
    the_parameters = parameter_defaults
    if kwargs:
        the_parameters = deepcopy(parameter_defaults)
        the_parameters.update(kwargs)
    quadratic_component = the_parameters['quadratic_component']
    new_objective = the_parameters['new_objective']
    if quadratic_component is not None:
        raise Exception('%s cannot solve QPs, try a different solver'%solver_name)
    #Faster to use these dicts than index lists
    index_to_metabolite = dict(zip(range(len(cobra_model.metabolites)),
                                   cobra_model.metabolites))
    index_to_reaction = dict(zip(range(len(cobra_model.reactions)),
                                 cobra_model.reactions))
    reaction_to_index = dict(zip(index_to_reaction.values(),
                                 index_to_reaction.keys()))


    lp = __solver_class()        # Create empty problem instance
    lp.name = 'cobra'     # Assign symbolic name to problem
    lp.rows.add(len(cobra_model.metabolites))
    lp.cols.add(len(cobra_model.reactions))
    linear_constraints = []
    for r in lp.rows:
        the_metabolite = index_to_metabolite[r.index]
        r.name = the_metabolite.id
        b = float(the_metabolite._bound)
        c = sense_dict[the_metabolite._constraint_sense]
        if c == 'E':
            r.bounds = b, b     # Set metabolite to steady state levels
        elif c == 'L':
            r.bounds = None, b
        elif c == 'G':
            r.bounds = b, None
        #Add in the linear constraints

        for the_reaction in the_metabolite._reaction:
            reaction_index = reaction_to_index[the_reaction]
            the_coefficient = the_reaction._metabolites[the_metabolite]
            linear_constraints.append((r.index, reaction_index,
                                       the_coefficient))
    #Need to assign lp.matrix after constructing the whole list
    lp.matrix = linear_constraints
    objective_coefficients = []

    for c in lp.cols:
        the_reaction = index_to_reaction[c.index]
        c.name = the_reaction.id           
        the_reaction = index_to_reaction[c.index]
        c.kind = variable_kind_dict[the_reaction.variable_kind]
        c.bounds = the_reaction.lower_bound, the_reaction.upper_bound
        objective_coefficients.append(float(the_reaction.objective_coefficient))
    #Add the new objective coefficients to the problem
    lp.obj[:] = objective_coefficients

    # make sure the objective sense is set in create_problem
    if "objective_sense" in the_parameters:
        set_parameter(lp, "objective_sense", the_parameters["objective_sense"])

    return(lp)


def update_problem(lp, cobra_model, **kwargs):
    """A performance tunable method for updating a model problem file

    lp: A gurobi problem object

    cobra_model: the cobra.Model corresponding to 'lp'

    """
    #When reusing the basis only assume that the objective coefficients or bounds can change
    #BUG with changing / unchanging the basis
    index_to_metabolite = dict(zip(range(len(cobra_model.metabolites)),
                                   cobra_model.metabolites))
    index_to_reaction = dict(zip(range(len(cobra_model.reactions)),
                                 cobra_model.reactions))
    reaction_to_index = dict(zip(index_to_reaction.values(),
                                 index_to_reaction.keys()))

    try:
        new_objective = kwargs['new_objective']
    except:
        new_objective = None
    if new_objective is not None:
        objective_coefficients = []
        for c in lp.cols:      # Iterate over all rows
            the_reaction = index_to_reaction[c.index]
            c.name = the_reaction.id
            c.bounds = the_reaction.lower_bound, the_reaction.upper_bound
            objective_coefficients.append(float(the_reaction.objective_coefficient))
            c.kind = variable_kind_dict[the_reaction.variable_kind]
        #Add the new objective coefficients to the problem
        lp.obj[:] = objective_coefficients
    else:
        for c in lp.cols:      # Iterate over all rows
            the_reaction = index_to_reaction[c.index]
            c.name = the_reaction.id
            c.bounds = the_reaction.lower_bound, the_reaction.upper_bound
            c.kind = variable_kind_dict[the_reaction.variable_kind]


def solve_problem(lp, **kwargs):
    """A performance tunable method for updating a model problem file

    lp: a pyGLPK 0.3 problem

    For pyGLPK it is necessary to provide the following parameters, if they
    are not provided then the default settings will be used: tolerance_optimality,
    tolerance_integer, lp_method, and objective_sense

    """
    if kwargs:
        the_parameters = kwargs
    else:
        the_parameters = {}
    #These need to be provided to solve_problem for pyGLPK because it's not possible
    #AFAIK to set these during problem creation.
    __function_parameters = ['tolerance_optimality', 'tolerance_integer', 'lp_method']
    for the_parameter in __function_parameters:
        if the_parameter not in the_parameters:
            the_parameters.update({the_parameter: parameter_defaults[the_parameter]})
           
    
    tolerance_optimality = the_parameters['tolerance_optimality']
    tolerance_integer = the_parameters['tolerance_integer']
    lp_method = the_parameters['lp_method']
    #[set_parameter(lp, parameter_mappings[k], v)
    # for k, v in kwargs.iteritems() if k in parameter_mappings]
    if "objective_sense" in the_parameters:
        set_parameter(lp, "objective_sense", the_parameters["objective_sense"])
    if lp.kind == int:
        #For MILPs, it is faster to solve LP then move to MILP
        lp.simplex(tol_bnd=tolerance_optimality,
                   tol_dj=tolerance_optimality, meth=lp_method)  
        lp.integer(tol_int=tolerance_integer)
    else:
        lp.simplex(tol_bnd=tolerance_optimality, tol_dj=tolerance_optimality,
                   meth=lp_method)
    status = get_status(lp)

    return status


def solve(cobra_model, **kwargs):
    """Smart interface to optimization solver functions that will convert
    the cobra_model to a solver object, set the parameters, and try multiple
    methods to get an optimal solution before returning the solver object and
    a cobra.Solution (which is attached to cobra_model.solution)

    cobra_model: a cobra.Model

    returns a dict: {'the_problem': solver specific object, 'the_solution':
    cobra.Solution for the optimization problem'}


    """
    #Start out with default parameters and then modify if
    #new onese are provided
    the_parameters = deepcopy(parameter_defaults)
    if kwargs:
        the_parameters.update(kwargs)
    #Update objectives if they are new.
    error_reporting = the_parameters['error_reporting']
    if 'new_objective' in the_parameters and \
           the_parameters['new_objective'] not in ['update problem', None]:
       update_objective(cobra_model, the_parameters['new_objective'])
    if 'the_problem' in the_parameters:
        the_problem = the_parameters['the_problem']
    else:
        the_problem = None
    if isinstance(the_problem, __solver_class):
        #Update the problem with the current cobra_model
        lp = the_problem
        update_problem(lp, cobra_model, **the_parameters)
    else:
        #Create a new problem
        lp = create_problem(cobra_model, **the_parameters)

    ###Try to solve the problem using other methods if the first method doesn't work
    lp_method = the_parameters['lp_method']
    the_methods = [1, 2, 3]
    if lp_method in the_methods:
        the_methods.remove(lp_method)
    #Start with the user specified method
    the_methods.insert(0, lp_method)
    for the_method in the_methods:
        the_parameters['lp_method'] = the_method
        try:
            status = solve_problem(lp, **the_parameters)
        except:
            status = 'failed'
        if status == 'optimal':
            break

    the_solution = format_solution(lp, cobra_model)
    if status != 'optimal' and error_reporting:
        print '%s failed: %s'%(solver_name, status)
    cobra_model.solution = the_solution
    solution = {'the_problem': lp, 'the_solution': the_solution}
    return solution
