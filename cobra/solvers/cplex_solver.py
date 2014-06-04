#cobra.solvers.cplex_solver
#Interface to ilog/cplex 12.4 python interface

from os import name as __name
from copy import deepcopy
from warnings import warn
###solver specific parameters
from .parameters import status_dict, variable_kind_dict, \
     sense_dict, parameter_mappings, parameter_defaults, \
     objective_senses, default_objective_sense

from ..core.Solution import Solution

from time import time
solver_name = 'cplex'
parameter_defaults = parameter_defaults[solver_name]
sense_dict = eval(sense_dict[solver_name])


parameter_mappings = parameter_mappings[solver_name]
objective_senses = objective_senses[solver_name]
from cplex import Cplex, SparsePair
class Problem(Cplex):
    def __init__(self):
        Cplex.__init__(self)
__solver_class = Problem
variable_kind_dict = eval(variable_kind_dict[solver_name])
status_dict = eval(status_dict[solver_name])
def get_status(lp):
    status = lp.solution.get_status_string().lower()
    if status in status_dict:
        status = status_dict[status]
    else:
        status = 'failed'
    return status

def get_objective_value(lp):
    return lp.solution.get_objective_value()

def format_solution(lp, cobra_model, **kwargs):
    status = get_status(lp)
    #TODO: It might be able to speed this up a little.
    if status in ('optimal', 'time_limit'):
        objective_value = lp.solution.get_objective_value()
        #This can be sped up a little
        x_dict = dict(zip(lp.variables.get_names(),
                     lp.solution.get_values()))
        x = lp.solution.get_values()
        #MIP's don't have duals
        if lp.get_problem_type() in (Cplex.problem_type.MIQP,
                                     Cplex.problem_type.MILP):

            y = y_dict = None
        else:
            y_dict = dict(zip(lp.linear_constraints.get_names(),
                              lp.solution.get_dual_values()))
            y = lp.solution.get_dual_values()
    else:
        x = y = x_dict = y_dict = objective_value = None

    the_solution = Solution(objective_value, x=x, x_dict=x_dict,
                            status=status, y=y, y_dict=y_dict)
    return the_solution    

def set_parameter(lp, parameter_name, parameter_value):
    """
    
    """
    if parameter_name == 'objective.set_sense':
        if parameter_value in objective_senses:
            parameter_value = eval(objective_senses[parameter_value])
    try:
        eval('lp.%s(%s)'%(parameter_name, repr(parameter_value)))
    except Exception, e:
        print "Couldn't set parameter %s: %s"%(parameter_name, repr(e))



def create_problem(cobra_model, quadratic_component=None, **kwargs):
    """Solver-specific method for constructing a solver problem from
    a cobra.Model.  This can be tuned for performance using kwargs


    """
    the_parameters = parameter_defaults
    if kwargs:
        the_parameters = deepcopy(parameter_defaults)
        the_parameters.update(kwargs)

    lp = __solver_class()
    if 'log_file' not in the_parameters:
        lp.set_results_stream(None)
        lp.set_warning_stream(None)
    [set_parameter(lp, parameter_mappings[k], v)
     for k, v in the_parameters.iteritems() if k in parameter_mappings]
    if 'relax_b' in the_parameters:
        warn('need to reimplement relax_b')
        relax_b = False
    else:
        relax_b = False

    #Using the new objects
    #NOTE: This might be slow
    objective_coefficients = []
    lower_bounds = []
    upper_bounds = []
    variable_names = []
    variable_kinds = []
    [(objective_coefficients.append(x.objective_coefficient),
      lower_bounds.append(x.lower_bound),
      upper_bounds.append(x.upper_bound),
      variable_names.append(x.id),
      variable_kinds.append(variable_kind_dict[x.variable_kind]))
     for x in cobra_model.reactions]
    #Cplex decides that the problem is a MIP if variable_kinds are supplied
    #even if there aren't any integers.
    if variable_kind_dict['integer'] in variable_kinds:
        lp.variables.add(obj=objective_coefficients,
                         lb=lower_bounds,
                         ub=upper_bounds,
                         names=variable_names,
                         types=variable_kinds)
    else:
        lp.variables.add(obj=objective_coefficients,
                         lb=lower_bounds,
                         ub=upper_bounds,
                         names=variable_names)

   ## if relax_b:
        ## range_values = zeros(len(cobra_model.metabolites))
        ## b_values = array([x._bound for x in cobra_model.metabolties])
        ## for the_nonzero in list(b_values.nonzero()[0]):
        ##     range_values[the_nonzero] = -relax_b

    constraint_sense = []
    constraint_names = []
    constraint_limits = []
    [(constraint_sense.append(x._constraint_sense),
      constraint_names.append(x.id),
      constraint_limits.append(x._bound))
     for x in cobra_model.metabolites]

    the_linear_expressions = []
    #NOTE: This won't work with metabolites that aren't in any reaction
    for the_metabolite in cobra_model.metabolites:
        variable_list = []
        coefficient_list = []
        for the_reaction in the_metabolite._reaction:
            variable_list.append(the_reaction.id)
            coefficient_list.append(the_reaction._metabolites[the_metabolite])
        the_linear_expressions.append(SparsePair(ind=variable_list,
                                                 val=coefficient_list))
    # Set objective to quadratic program
    if quadratic_component is not None:
        set_quadratic_objective(lp, quadratic_component)

    if relax_b:
        lp.linear_constraints.add(lin_expr=the_linear_expressions,
                                  rhs=constraint_limits,
                                  range_values=list(range_values),
                                  senses=constraint_sense,
                                  names=constraint_names)

    else:
        lp.linear_constraints.add(lin_expr=the_linear_expressions,
                                  rhs=constraint_limits,
                                  senses=constraint_sense,
                                  names=constraint_names)

    #Set the problem type as cplex doesn't appear to do this correctly
    problem_type = Cplex.problem_type.LP
    if Cplex.variables.type.integer in variable_kinds:
        if quadratic_component is not None:
            problem_type = Cplex.problem_type.MIQP
        else:
            problem_type = Cplex.problem_type.MILP
    elif quadratic_component is not None:
        problem_type = Cplex.problem_type.QP
    lp.set_problem_type(problem_type)
    return(lp)


def set_quadratic_objective(lp, quadratic_objective):
    if not hasattr(quadratic_objective, 'todok'):
        raise Exception('quadratic component must have method todok')

    # Reset the quadratic coefficient if it exists
    if lp.objective.get_num_quadratic_nonzeros() > 0:
        lp.objective.set_quadratic((0.,) * lp.variables.get_num())

    # cplex divides by 2 for some reason
    quadratic_component_scaled = quadratic_objective.todok() * 2

    lp.parameters.emphasis.numerical.set(1)
    for k, v in quadratic_component_scaled.items():
        lp.objective.set_quadratic_coefficients(int(k[0]), int(k[1]), v)

def change_variable_bounds(lp, index, lower_bound, upper_bound):
    lp.variables.set_lower_bounds(index, lower_bound)
    lp.variables.set_upper_bounds(index, upper_bound)

def change_variable_objective(lp, index, objective):
    lp.objective.set_linear(index, objective)


def change_coefficient(lp, met_index, rxn_index, value):
    lp.linear_constraints.set_coefficients(met_index, rxn_index, value)

def update_problem(lp, cobra_model, **kwargs):
    """A performance tunable method for updating a model problem file

    lp: A cplex problem object

    cobra_model: the cobra.Model corresponding to 'lp'

    """
    #When reusing the basis only assume that the objective coefficients or bounds can change
    #BUG with changing / unchanging the basis

    try:
        update_problem_reaction_bounds = kwargs['update_problem_reaction_bounds']
    except:
        update_problem_reaction_bounds = True
    try:
        quadratic_component = kwargs['quadratic_component']
        if quadratic_component is not None:
            warn("update_problem does not yet take quadratic_component as a parameter")
    except:
        quadratic_component = None

    if new_objective is not None:
        lp.objective.set_linear([(x.id, float(x.objective_coefficient))
                                 for x in cobra_model.reactions])
    if update_problem_reaction_bounds:
        lp.variables.set_upper_bounds([(x.id, float(x.upper_bound))
                                        for x in cobra_model.reactions])
        lp.variables.set_lower_bounds([(x.id, float(x.lower_bound))
                                        for x in cobra_model.reactions])



###
def solve_problem(lp, **kwargs):
    """A performance tunable method for solving a problem

    """
    #Update parameter settings if provided
    if kwargs:
        [set_parameter(lp, parameter_mappings[k], v)
         for k, v in kwargs.iteritems() if k in parameter_mappings]
    try:
        the_problem = kwargs['the_problem']
    except:
        the_problem = False
    if isinstance(the_problem, __solver_class):
        try:
            the_basis = the_problem.solution.basis.get_basis()
            lp.start.set_basis(the_basis[0],the_basis[1])
            lp.parameters.preprocessing.presolve.set(0)
        except:
            warn("cplex_java isn't yet configured to reuse the basis")
    
    lp.solve()
    #If the solver takes more than 0.1 s with a hot start it is likely stuck
    status = get_status(lp)
    return status

    
def solve(cobra_model, **kwargs):
    """

    """
    #Start out with default parameters and then modify if
    #new onese are provided
    the_parameters = deepcopy(parameter_defaults)
    if kwargs:
        the_parameters.update(kwargs)
    for i in ["new_objective", "update_problem", "the_problem"]:
        if i in the_parameters:
            raise Exception("Option %s removed" % i)
    if 'error_reporting' in the_parameters:
        warn("error_reporting deprecated")


    lp = create_problem(cobra_model, **the_parameters)

    ###Try to solve the problem using other methods if the first method doesn't work
    try:
        lp_method = the_parameters['lp_method']
    except:
        lp_method = 1
    the_methods = [1, 2, 3, 4, 5, 6]
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
    #if status != 'optimal':
    #    print '%s failed: %s'%(solver_name, status)
    #cobra_model.solution = the_solution
    #solution = {'the_problem': lp, 'the_solution': the_solution}
    return the_solution
