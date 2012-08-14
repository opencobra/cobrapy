# because both this module and the library are named cplex
from __future__ import absolute_import

from cplex import Cplex, SparsePair
from ..core.Solution import Solution
from scipy.sparse import dok_matrix
from numpy import array
from sys import stdout, stderr

# sense_dict = {}
variable_kind_dict = {
    'continuous': Cplex.variables.type.continuous,
    'integer': Cplex.variables.type.integer}
parameter_dict = {## todo - find a way to use
    "threads": "threads",
    "timelimit": "timelimit",
    "MIPGapAbs": "mip.tolerances.absmipgap",
    "MIPGap": "mip.tolerances.mipgap"
}
default_parameters = { ## todo find a way to sue
}

status_dict = {
    'MIP_infeasible': 'infeasible',
    'MIP_optimal': 'optimal',
    'MIP_optimal_tolerance': 'optimal',
    'MIP_unbounded': 'unbounded',
    'infeasible': 'infeasible',
    'optimal': 'optimal',
    'optimal_tolerance': 'optimal',
    'integer optimal solution': 'optimal',
    'unbounded': 'unbounded',
    'integer optimal, tolerance': 'optimal',
    'time limit exceeded': 'time_limit'}

def create_problem(cobra_model, objective_sense="maximize",
        quadratic_component=None):
    """set up the CPLEX solver object
    
    Parameters
    cobra_model: :class:`~cobra.core.Model`
    objective_sense: either "maximize" or "minimize"
    quadratic_component: matrix representing the quadratic objective function
    """
    # TODO - implement relax_b
    lp = Cplex()
    lp.set_results_stream(None)
    lp.set_warning_stream(None)
    objective_coefficients = []
    lower_bounds = []
    upper_bounds = []
    variable_names = []
    variable_kinds = []
    for x in cobra_model.reactions:
        objective_coefficients.append(float(x.objective_coefficient))
        lower_bounds.append(float(x.lower_bound))
        upper_bounds.append(float(x.upper_bound))
        variable_names.append(str(x.id))
        variable_kinds.append(variable_kind_dict[x.variable_kind])
    #Cplex decides that the problem is a MIP if variable_kinds are supplied
    #even if there aren't any integers.
    # This can probably be skipped, because we manually set the problem kind later
    if Cplex.variables.type.integer in variable_kinds:
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
        constraint_sense.append(str(x._constraint_sense))
        constraint_names.append(str(x.id))
        constraint_limits.append(float(x._bound))
    
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
    if quadratic_component is not None:
        quadratic_component_scaled = dok_matrix(quadratic_component)
        lp.parameters.emphasis.numerical.set(1)
        for k, v in quadratic_component_scaled.items():
            lp.objective.set_quadratic_coefficients(int(k[0]), int(k[1]), v)

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
    if objective_sense == 'maximize':
        lp.objective.set_sense(lp.objective.sense.maximize)
    elif objective_sense == 'minimize':
        lp.objective.set_sense(lp.objective.sense.minimize)
    else:
        raise ValueError("objective_sense")
    return lp

def _set_parameter(lp, parameter_name, parameter_value):
    parameter = lp.parameters
    for attribute in parameter_name.split("."):
        parameter = getattr(parameter, attribute)
    parameter.set(parameter_value)


def solve_problem(lp, **parameters):
    """
    Parameters
    threads: number of threads to use
    verbose: Boolean"""
    if "verbose" in parameters:
        if parameters["verbose"]:
            lp.set_results_stream(stdout)
            lp.set_warning_stream(stdout)
        parameters.pop("verbose")
    for key, value in parameters.iteritems():
        if key in parameter_dict:
            _set_parameter(lp, parameter_dict[key], value)
        else:
            raise Exception("parameter %s not understood" % key)
    # if "threads" in parameters:
        # lp.parameters.threads.set(parameters["threads"])
    # if "timelimit" in parameters:
        # lp.parameters.set(parameters["timelimit"])
    lp.solve()
    status = status_dict[lp.solution.get_status_string()]
    solution = Solution(status)
    solution.status = status
    if status == "optimal" or status == "time_limit":
        solution.objective_value = lp.solution.get_objective_value()
        #This can be sped up a little
        solution.x_dict = dict(zip(lp.variables.get_names(),
                     lp.solution.get_values()))
        # solution.x = array(lp.solution.get_values())
        solution.x = lp.solution.get_values()
        #MIP's don't have duals
        if lp.get_problem_type() in (Cplex.problem_type.MIQP,
                                     Cplex.problem_type.MILP):
            solution.y = solution.y_dict = None
        else:
            solution.y_dict = dict(zip(lp.linear_constraints.get_names(),
                              lp.solution.get_dual_values()))
            y = array(lp.solution.get_dual_values())
            solution.y = y.reshape(y.shape[0],1)
    return solution


def solve(cobra_model, objective_sense="maximize", **kwargs):
    return solve_problem(create_problem(cobra_model, objective_sense=objective_sense), **kwargs)

def _optimize_cplex(cobra_model, new_objective=None, objective_sense='maximize',
                   min_norm=0, the_problem=None, 
                   tolerance_optimality=1e-6, tolerance_feasibility=1e-6, tolerance_integer=1e-9,
                   tolerance_barrier=1e-8,error_reporting=None, 
                   print_solver_time=False, lp_method=1, lp_parallel=0, copy_problem=False,
                   relax_b=None, quadratic_component=None, reuse_basis=True,
                   update_problem_reaction_bounds=True):
   

    if not error_reporting:
        lp.set_results_stream(None)
        lp.set_warning_stream(None)
    if  print_solver_time:
        start_time = time()
    if not isinstance(the_problem, Cplex):
        #TODO: set tolerance
        lp.solve()
        # Solve this LP with the simplex method.  Takes about 0.2 s without hot start
        lp.status = lp.solution.status[lp.solution.get_status()]
        if lp.status in status_dict:
            status = status_dict[lp.status]
        else:
            status = 'failed'
    else:
        if isinstance(the_problem, Cplex) and reuse_basis:
            try:
                the_basis = the_problem.solution.basis.get_basis()
                lp.start.set_basis(the_basis[0],the_basis[1])
                #TODO: Determine whether the primal or dual works best for the
                #problem of interest.  For the ME matrix the primal appears to
                #work best
                lp_method = 1
                lp.parameters.preprocessing.presolve.set(0)
                lp.parameters.lpmethod.set(lp_method) 
            except:
                print 'no basis in the_problem'
        #TODO: set tolerance and time limit
        #lp.parameters.timelimit.set()
        lp.solve()
        #If the solver takes more than 0.1 s with a hot start it is likely stuck
        lp.status = lp.solution.status[lp.solution.get_status()]
        if lp.status in status_dict:
            status = status_dict[lp.status]
        else:
            status = 'failed'
        if status != 'optimal':
            #Cycle through the different solver options, if a solution is not found
            for lp_method in (1, 2, 3, 4, 5, 6):
                lp = optimize_cplex(cobra_model, new_objective=new_objective,
                                    objective_sense=objective_sense,
                                    min_norm=min_norm, the_problem=None,
                                    print_solver_time=print_solver_time,
                                    tolerance_optimality=tolerance_optimality,
                                    tolerance_feasibility=tolerance_feasibility,
                                    lp_method=lp_method,
                                    quadratic_component=quadratic_component)['the_problem']
                lp.status = lp.solution.status[lp.solution.get_status()]
                if lp.status in status_dict:
                    status = status_dict[lp.status]
                else:
                    status = 'failed'
                if status == 'optimal':
                    break
    if error_reporting == 'time':
        print 'solver time: ' + repr(time()-start_time) + ' with method ' + repr(lp_method)
        start_time = time()

    if print_solver_time:
        print 'cplex time: %f'%(time() - start_time)
    x = []
    x_dict = {}
    #TODO: It might be able to speed this up a little.
    if status == 'optimal':
        objective_value = lp.solution.get_objective_value()
        #This can be sped up a little
        x_dict = dict(zip(lp.variables.get_names(),
                     lp.solution.get_values()))
        x = array(lp.solution.get_values())
        x = x.reshape(x.shape[0],1)
        #MIP's don't have duals
        if lp.get_problem_type() in (Cplex.problem_type.MIQP,
                                     Cplex.problem_type.MILP):

            y = y_dict = None
        else:
            y_dict = dict(zip(lp.linear_constraints.get_names(),
                              lp.solution.get_dual_values()))
            y = array(lp.solution.get_dual_values())
            y = y.reshape(y.shape[0],1)
    else:
        x = y = x_dict = y_dict = objective_value = None
        if error_reporting:
            print 'cplex failed: %s'%lp.status

    the_solution = Solution(objective_value, x=x, x_dict=x_dict,
                            status=status, y=y, y_dict=y_dict)
    solution = {'the_problem': lp, 'the_solution': the_solution}
    return solution    