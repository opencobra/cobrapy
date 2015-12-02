# Interface to ilog/cplex 12.4 python interface

from copy import deepcopy
from warnings import warn
import sys

from cplex import Cplex, SparsePair
from cplex.exceptions import CplexError

from ..core.Solution import Solution
from six.moves import zip
from six import string_types, iteritems

try:
    from sympy import Basic, Number
except:
    class Basic:
        pass

def _float(value):
    if isinstance(value, Basic) and not isinstance(value, Number):
        return 0.
    else:
        return float(value)

solver_name = 'cplex'
_SUPPORTS_MILP = True

# solver specific parameters
parameter_defaults = {'objective_sense': 'maximize',
                      'tolerance_optimality': 1e-6,
                      'tolerance_feasibility': 1e-6,
                      'tolerance_integer': 1e-9,
                      'lp_method': 1,
                      'tolerance_barrier': 1e-8,
                      'verbose': False,
                      'qpmethod': 1}
parameter_mappings = {'lp_method': 'lpmethod',
                      'lp_parallel': 'threads',
                      'threads': 'threads',
                      'objective_sense': 'objective_sense',
                      'time_limit': 'timelimit',
                      'iteration_limit': 'simplex.limits.iterations',
                      'tolerance_barrier': 'barrier.convergetol',
                      'tolerance_feasibility': 'simplex.tolerances.feasibility',
                      'tolerance_markowitz': 'simplex.tolerances.markowitz',
                      'tolerance_optimality': 'simplex.tolerances.optimality',
                      'tolerance_integer': 'mip.tolerances.integrality',
                      'MIP_gap_abs': 'mip.tolerances.absmipgap',
                      'MIP_gap': 'mip.tolerances.mipgap'}
variable_kind_dict = {'continuous': Cplex.variables.type.continuous,
                      'integer': Cplex.variables.type.integer}
status_dict = {'MIP_infeasible': 'infeasible',
               'integer optimal solution': 'optimal',
               'MIP_optimal': 'optimal',
               'MIP_optimal_tolerance': 'optimal',
               'MIP_unbounded':  'unbounded',
               'infeasible': 'infeasible',
               'integer infeasible': 'infeasible',
               'optimal': 'optimal',
               'optimal_tolerance': 'optimal',
               'unbounded': 'unbounded',
               'integer optimal, tolerance': 'optimal',
               'time limit exceeded': 'time_limit'}


def get_status(lp):
    status = lp.solution.get_status_string().lower()
    return status_dict.get(status, status)


def get_objective_value(lp):
    return lp.solution.get_objective_value()


def format_solution(lp, cobra_model, **kwargs):
    status = get_status(lp)
    if status in ('optimal', 'time_limit', 'non-optimal'):
        objective_value = lp.solution.get_objective_value()
        x_dict = dict(zip(lp.variables.get_names(),
                          lp.solution.get_values()))
        x = lp.solution.get_values()
        # MIP's don't have duals
        if lp.get_problem_type() in (Cplex.problem_type.MIQP,
                                     Cplex.problem_type.MILP):

            y = y_dict = None
        else:
            y_dict = dict(zip(lp.linear_constraints.get_names(),
                              lp.solution.get_dual_values()))
            y = lp.solution.get_dual_values()
    else:
        x = y = x_dict = y_dict = objective_value = None

    return Solution(objective_value, x=x, x_dict=x_dict, status=status,
                    y=y, y_dict=y_dict)


def set_parameter(lp, parameter_name, parameter_value):
    if parameter_name == 'objective_sense':
        parameter_value = getattr(lp.objective.sense, parameter_value)
        lp.objective.set_sense(parameter_value)
        return
    elif parameter_name == 'the_problem':
        warn('option the_problem removed')
        return
    elif parameter_name == 'verbose':
        if parameter_value:
            lp.set_log_stream(sys.stdout)
            lp.set_results_stream(sys.stdout)
            lp.set_warning_stream(sys.stderr)
            # If the value passed in is True, it shold be 1. MIP display can
            # be as high as 5, but the others only go up to 2.
            value = int(parameter_value)
            set_parameter(lp, 'mip.display', value)
            set_parameter(lp, 'simplex.display', min(value, 2))
            set_parameter(lp, 'barrier.display', min(value, 2))
        else:
            lp.set_log_stream(None)
            lp.set_results_stream(None)
            lp.set_warning_stream(None)
            set_parameter(lp, 'mip.display', 0)
            set_parameter(lp, 'simplex.display', 0)
            set_parameter(lp, 'barrier.display', 0)
        return
    try:
        cplex_name = parameter_mappings.get(parameter_name, parameter_name)
        cplex_value = parameter_value
        # This will iteratively get parameters. For example
        # "simplex.tolerances.feasibility" will evaluate to
        # lp.parameters.simplex.tolerances.feasibility.
        param = lp.parameters
        for i in cplex_name.split("."):
            param = getattr(param, i)
        # For example, this next part will allow setting the parameter
        # lpmethod to "auto", "primal", "dual" or any other string in
        # parameters.lp.method.values
        if isinstance(cplex_value, string_types) and \
                hasattr(param.values, cplex_value):
            cplex_value = getattr(param.values, cplex_value)
        param.set(cplex_value)
    except (CplexError, AttributeError) as e:
        raise ValueError("Failed to set %s to %s: %s" %
                         (parameter_name, str(parameter_value), repr(e)))


def create_problem(cobra_model, quadratic_component=None, **kwargs):
    """Solver-specific method for constructing a solver problem from
    a cobra.Model.  This can be tuned for performance using kwargs


    """
    # Process parameter defaults
    the_parameters = parameter_defaults
    if kwargs:
        the_parameters = parameter_defaults.copy()
        the_parameters.update(kwargs)
    if 'relax_b' in the_parameters:
        relax_b = the_parameters.pop("relax_b")
        warn('need to reimplement relax_b')
        relax_b = False
    else:
        relax_b = False

    # Begin problem creation
    lp = Cplex()
    for k, v in iteritems(the_parameters):
        set_parameter(lp, k, v)
    objective_coefficients = [float(x.objective_coefficient)
                              for x in cobra_model.reactions]
    lower_bounds = [_float(x.lower_bound) for x in cobra_model.reactions]
    upper_bounds = [_float(x.upper_bound) for x in cobra_model.reactions]
    variable_names = cobra_model.reactions.list_attr("id")
    variable_kinds = [variable_kind_dict[x.variable_kind] for x
                      in cobra_model.reactions]
    # Cplex decides that the problem is a MIP if variable_kinds are supplied
    # even if there aren't any integers.
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

    constraint_sense = []
    constraint_names = []
    constraint_limits = []

    for x in cobra_model.metabolites:
        constraint_sense.append(x._constraint_sense)
        constraint_names.append(x.id)
        constraint_limits.append(float(x._bound))

    the_linear_expressions = []
    # NOTE: This won't work with metabolites that aren't in any reaction
    for the_metabolite in cobra_model.metabolites:
        variable_list = []
        coefficient_list = []
        for the_reaction in the_metabolite._reaction:
            variable_list.append(the_reaction.id)
            coefficient_list.append(_float(the_reaction._metabolites[the_metabolite]))
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

    # Set the problem type as cplex doesn't appear to do this correctly
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
    # ensure the matrix is properly read in
    nnz = quadratic_objective.nnz
    if lp.parameters.read.qpnonzeros.get() < nnz:
        lp.parameters.read.qpnonzeros.set(nnz + 1)
    # Reset the quadratic coefficient if it exists
    if lp.objective.get_num_quadratic_nonzeros() > 0:
        lp.objective.set_quadratic((0.,) * lp.variables.get_num())

    quadratic_component_scaled = quadratic_objective.todok()

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


def update_problem(lp, cobra_model, new_objective=None, **kwargs):
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


def solve_problem(lp, **kwargs):
    # Update parameter settings if provided
    for k, v in iteritems(kwargs):
        set_parameter(lp, k, v)
    lp.solve()
    # If the solver takes more than 0.1 s with a hot start it is likely stuck
    return get_status(lp)


def solve(cobra_model, **kwargs):
    """

    """
    # Start out with default parameters and then modify if
    # new onese are provided
    for i in ["new_objective", "update_problem", "the_problem"]:
        if i in kwargs:
            raise Exception("Option %s removed" % i)
    if 'error_reporting' in kwargs:
        kwargs.pop('error_erporting')
        warn("error_reporting deprecated")

    # Create problem will get parameter defaults
    lp = create_problem(cobra_model, **kwargs)

    # Try to solve the problem using other methods if first method doesn't work
    try:
        lp_method = the_parameters['lp_method']
    except:
        lp_method = 1
    the_methods = [1, 2, 3, 4, 5, 6]
    if lp_method in the_methods:
        the_methods.remove(lp_method)
    # Start with the user specified method
    the_methods.insert(0, lp_method)
    for the_method in the_methods:
        try:
            status = solve_problem(lp, lp_method=the_method)
        except:
            status = 'failed'
        if status == 'optimal':
            break
    return format_solution(lp, cobra_model)
