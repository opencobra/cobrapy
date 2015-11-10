# PLEASE NOTE THAT JYTHON SUPPORT (and this jython-only-solver) is deprecated
#Interface to ilog/cplex 12.4 python / jython interfaces
#QPs are not yet supported under jython
from __future__ import print_function
from os import name as __name
from copy import deepcopy
from warnings import warn
###solver specific parameters
from .parameters import status_dict, variable_kind_dict, \
     sense_dict, parameter_mappings, parameter_defaults, \
     objective_senses, default_objective_sense

from ..core.Solution import Solution
from time import time
from six import iteritems
solver_name = 'cplex'
parameter_defaults = parameter_defaults[solver_name]
sense_dict = eval(sense_dict[solver_name])

from ilog.cplex import IloCplex
from ilog.cplex.IloCplex import DoubleParam, IntParam, StringParam
from ilog.concert import IloNumVarType, IloObjectiveSense 
#__solver_class = IloCplex
status_dict = eval(status_dict[solver_name])
class Problem(IloCplex):
    def __init__(self):
        IloCplex.__init__(self)
        self._lp_matrix = self.addLPMatrix()
        self.objective_value = None
        self._objective_sense = 'maximize'
    def add_linear_expression(self, linear_expression, metabolite):
        b = metabolite._bound
        c = metabolite._constraint_sense
        the_id = metabolite.id
        if c == 'E':
            p = self.addEq(linear_expression, b, the_id)
        elif c == 'L':
            p = self.addLe(linear_expression, b, the_id)
        elif c == 'G':
            p = self.addGe(linear_expression, b, the_id)
        else:
            raise Exception("Constraint sense '%s' for metabolite %s is not valid"%(c,
                                                                                    the_id))
        return(p)

__solver_class = Problem
parameter_mappings = parameter_mappings['%s_%s'%(solver_name,
                                                 __name)]
variable_kind_dict = eval(variable_kind_dict['%s_%s'%(solver_name,
                                                      __name)])
objective_senses = objective_senses['%s_%s'%(solver_name,
                                             __name)]
## from jarray import array as j_array
## def array(x, variable_type='d'):
##     return j_array(x, variable_type)


def get_status(lp):
    status = repr(lp.status).lower()
    if status in status_dict:
        status = status_dict[status]
    else:
        status = 'failed'
    return status

def set_parameter(lp, parameter_name, parameter_value):
    if parameter_name == 'objective_sense':
        objective = lp.getObjective()
        if objective is not None:
            objective.setSense(eval(objective_senses[parameter_value]))

    else:
        if hasattr(DoubleParam, parameter_name):
            parameter_type = 'DoubleParam'
        elif hasattr(IntParam, parameter_name):
            parameter_type = 'IntParam'
        elif hasattr(StringParam, parameter_name):
            parameter_type = 'StringParam'
        else:
            raise Exception("%s is not a DoubleParam, IntParam, or StringParam"%parameter_name)
        lp.setParam(eval('%s.%s'%(parameter_type, parameter_name)),
                    parameter_value)

def get_objective_value(lp):
    return lp.getObjValue()

def format_solution(lp, cobra_model, **kwargs):
    """
    TODO
    """
    status = get_status(lp)
    try:
        x = lp.getValues(lp.variables)
        x_dict = dict(zip(cobra_model.reactions, x))
        objective_value = lp.getObjValue()
    except:
        x = x_dict = objective_value = None
        #print status

    try:
        y = lp.getDuals(lp.variables)
        y_dict = dict(zip(cobra_model.metabolites, y))
    except:
        y = y_dict = None
    return(Solution(objective_value, x=x, x_dict=x_dict, y=y,
                    y_dict=y_dict, status=status))

def create_problem(cobra_model,  **kwargs):
    """Solver-specific method for constructing a solver problem from
    a cobra.Model.  This can be tuned for performance using kwargs

    TODO: This will need to be specific for python / jython
    """
    the_parameters = parameter_defaults
    if kwargs:
        the_parameters = deepcopy(parameter_defaults)
        the_parameters.update(kwargs)

    lp = Problem()

    if 'log_file' not in the_parameters:
        lp.setWarning(None)
        lp.setOut(None)
    [set_parameter(lp, parameter_mappings[k], v)
     for k, v in iteritems(the_parameters) if k in parameter_mappings]
    quadratic_component = the_parameters['quadratic_component']
    new_objective = the_parameters['new_objective']
    error_reporting = the_parameters['error_reporting']
    lp._objective_sense  = the_parameters['objective_sense']
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

    #Only add the variable types if one or more variables is an integer, just
    #in case the java interface has the same bug as the python interface where
    #the problem type switches to integer if variable types are supplied even
    #if all are continuous
    if variable_kind_dict['integer'] in variable_kinds:
        lp.variables = lp.numVarArray(len(cobra_model.reactions), lower_bounds,
                                      upper_bounds, variable_kinds, variable_names)
    else:
        lp.variables = lp.numVarArray(len(cobra_model.reactions), lower_bounds, upper_bounds,
                                      variable_names)
    
    lp.variable_dict = dict(zip(cobra_model.reactions, lp.variables))
    if lp._objective_sense == 'maximize':
        __lp_add_objective = lp.addMaximize
    else:
        __lp_add_objective = lp.addMinimize

    __lp_add_objective(lp.scalProd(lp.variables, objective_coefficients))
    

    
    lp.constraints = []
    lp.constraint_dict = {}
    for the_metabolite in cobra_model.metabolites:
        linear_expression = lp.sum([lp.prod(k._metabolites[the_metabolite],
                                            lp.variable_dict[k])
                                    for k in the_metabolite._reaction])
        expression_pointer = lp.add_linear_expression(linear_expression, the_metabolite)
        lp.constraints.append(expression_pointer)
        lp.constraint_dict[the_metabolite] = expression_pointer
    
    if quadratic_component is not None:
        raise Exception("cplex through java isn't configured for QPs yet")
        if not hasattr(quadratic_component, 'todok'):
            raise Exception('quadratic component must have method todok')
        quadratic_component_scaled = quadratic_component.todok()

        lp.parameters.emphasis.numerical.set(1)
        for k, v in quadratic_component_scaled.items():
            lp.objective.set_quadratic_coefficients(int(k[0]), int(k[1]), v)

    ## #Set the problem type as cplex doesn't appear to do this correctly
    ## problem_type = Cplex.problem_type.LP
    ## if Cplex.variables.type.integer in variable_kinds:
    ##     if quadratic_component is not None:
    ##         problem_type = Cplex.problem_type.MIQP
    ##     else:
    ##         problem_type = Cplex.problem_type.MILP
    ## elif quadratic_component is not None:
    ##     problem_type = Cplex.problem_type.QP
    ## lp.set_problem_type(problem_type)
    return(lp)


def update_problem(lp, cobra_model, **kwargs):
    """A performance tunable method for updating a model problem file

    lp: A gurobi problem object

    cobra_model: the cobra.Model corresponding to 'lp'

    """
    #When reusing the basis only assume that the objective coefficients or bounds can change
    #BUG with changing / unchanging the basis
    try:
        new_objective = kwargs['new_objective']
    except:
        new_objective = None
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

    the_objective = lp.getObjective()
    for the_variable, the_reaction in zip(lp.variables, cobra_model.reactions):
        the_variable.setUB(float(the_reaction.upper_bound))
        the_variable.setLB(float(the_reaction.lower_bound))
        the_objective.setLinearCoef(the_variable, the_reaction.objective_coefficient)
            


###
def solve_problem(lp, **kwargs):
    """A performance tunable method for solving a problem

    """
    #Update parameter settings if provided
    if kwargs:
        [set_parameter(lp, parameter_mappings[k], v)
         for k, v in iteritems(kwargs) if k in parameter_mappings]
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
    #Update objectives if they are new.
    if 'new_objective' in the_parameters:
        raise ValueError("new_objective option removed")

    if 'the_problem' in the_parameters:
        the_problem = the_parameters['the_problem']
    else:
        the_problem = None
    if 'error_reporting' in the_parameters:
        error_reporting = the_parameters['error_reporting']
    else:
        error_reporting = False
    if isinstance(the_problem, __solver_class):
        #Update the problem with the current cobra_model
        lp = the_problem
        update_problem(lp, cobra_model, **the_parameters)
    else:
        #Create a new problem
        lp = create_problem(cobra_model, **the_parameters)
    #Deprecated way for returning a solver problem created from a cobra_model
    #without performing optimization
    if the_problem == 'setup':
        return lp

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
    if status != 'optimal' and error_reporting:
        print('{:s} failed: {:s}'.format(solver_name, status))
    cobra_model.solution = the_solution
    solution = {'the_problem': lp, 'the_solution': the_solution}
    return solution
