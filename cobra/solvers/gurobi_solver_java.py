# PLEASE NOTE THAT JYTHON SUPPORT (and this jython-only-solver) is deprecated
#Interface to the gurobi 5.0.1 python and java solvers
#QPs are not yet supported on java
from __future__ import print_function
from warnings import warn
from os import name as __name
from copy import deepcopy
from six import iteritems
###solver specific parameters
from .parameters import status_dict, variable_kind_dict, \
     sense_dict, parameter_mappings, parameter_defaults, \
     objective_senses, default_objective_sense

from ..core.Solution import Solution
from time import time
solver_name = 'gurobi'
objective_senses = objective_senses[solver_name]
parameter_mappings = parameter_mappings[solver_name]
parameter_defaults = parameter_defaults[solver_name]
#Functions that are different for java implementation of a solver
## from jarray import array as j_array
## def array(x, variable_type='d'):
##     return j_array(x, variable_type)

from gurobi import GRB
variable_kind_dict = eval(variable_kind_dict[solver_name])
status_dict = eval(status_dict[solver_name])

from gurobi import GRBModel, GRBEnv
from gurobi import GRBLinExpr
from gurobi import GRBQuadExpr as QuadExpr
__solver_class = GRBModel
#TODO: Create a pythonesqe class similar to in glpk_solver
def Model(name=''):
    grb_environment = GRBEnv(name)
    tmp_model = GRBModel(grb_environment)
    return tmp_model
def LinExpr(coefficients, variables):
    coefficients, variables = map(list, [coefficients, variables])
    tmp_expression = GRBLinExpr()
    tmp_expression.addTerms(coefficients, variables)
    return tmp_expression

def get_status(lp):
    status = lp.get(GRB.IntAttr.Status)
    if status in status_dict:
        status = status_dict[status]
    else:
        status = 'failed'
    return status

def set_parameter(lp, parameter_name, parameter_value):
    """Sets model parameters and attributes.
    
    """
    grb_environment = lp.getEnv()
    try:
        if hasattr(GRB.DoubleParam, parameter_name):
            grb_environment.set(eval('GRB.DoubleParam.%s'%parameter_name),
                                     parameter_value)
        elif hasattr(GRB.IntParam, parameter_name):
            grb_environment.set(eval('GRB.IntParam.%s'%parameter_name),
                                     parameter_value)
        elif hasattr(GRB.StringParam, parameter_name):
            grb_environment.set(eval('GRB.StringParam.%s'%parameter_name),
                                parameter_value)
        elif hasattr(GRB.IntAttr, parameter_name):
            if parameter_name == 'ModelSense':
                parameter_value = objective_senses[parameter_value]
            lp.set(eval('GRB.IntAttr.%s'%parameter_name),
                                parameter_value)
        else:
            warn("%s is not a DoubleParam, IntParam, StringParam, IntAttr"%parameter_name)
            ## raise Exception("%s is not a DoubleParam, IntParam, StringParam, IntAttr"%parameter_name)
    except Exception as e:
        warn("%s %s didn't work %s"%(parameter_name, parameter_value, e))

def get_objective_value(lp):
    return lp.get(GRB.DoubleAttr.ObjVal)

def format_solution(lp, cobra_model, **kwargs):
    """
    """
    status = get_status(lp)
    if status not in ('optimal', 'time_limit'):
        the_solution = Solution(None, status=status)
    else:
        x_dict = dict(((v.get(GRB.StringAttr.VarName),
                        v.get(GRB.DoubleAttr.X))
                       for v in lp.getVars()))
        x = [x_dict[v.id] for v in cobra_model.reactions]
        objective_value = lp.get(GRB.DoubleAttr.ObjVal)
        if lp.get(GRB.IntAttr.IsMIP) != 0:
            y = y_dict = None #MIP's don't have duals
        else:
            y_dict = dict(((c.get(GRB.StringAttr.ConstrName), c.get(GRB.DoubleAttr.Pi))
                          for c in lp.getConstrs()))
            y = list([y_dict[v.id] for v in cobra_model.metabolites])
        the_solution = Solution(objective_value, x=x, x_dict=x_dict, y=y,
                                y_dict=y_dict, status=status)
    return(the_solution)

def update_problem(lp, cobra_model, **kwargs):
    """A performance tunable method for updating a model problem file

    lp: A gurobi problem object

    cobra_model: the cobra.Model corresponding to 'lp'

    """
    #When reusing the basis only assume that the objective coefficients or bounds can change
    try:
        quadratic_component = kwargs['quadratic_component']
        if quadratic_component is not None:
            warn("update_problem does not yet take quadratic_component as a parameter")
    except:
        quadratic_component = None

    if 'reuse_basis' in kwargs and not kwargs['reuse_basis']:
        lp.reset()
    for the_variable, the_reaction in zip(lp.getVars(),
                                          cobra_model.reactions):
        the_variable.set(GRB.DoubleAttr.LB, float(the_reaction.lower_bound))
        the_variable.set(GRB.DoubleAttr.UB, float(the_reaction.upper_bound))
        the_variable.set(GRB.DoubleAttr.Obj, float(the_reaction.objective_coefficient))




###
sense_dict = eval(sense_dict[solver_name])
def create_problem(cobra_model,  **kwargs):
    """Solver-specific method for constructing a solver problem from
    a cobra.Model.  This can be tuned for performance using kwargs


    """
    lp = Model("")
    #Silence the solver
    set_parameter(lp, 'OutputFlag', 0)

    the_parameters = parameter_defaults
    if kwargs:
        the_parameters = deepcopy(parameter_defaults)
        the_parameters.update(kwargs)

    [set_parameter(lp, parameter_mappings[k], v)
         for k, v in iteritems(the_parameters) if k in parameter_mappings]
    quadratic_component = the_parameters['quadratic_component']
    objective_sense = objective_senses[the_parameters['objective_sense']]


    # Create variables
    #TODO:  Speed this up
    variable_list = [lp.addVar(float(x.lower_bound),
                               float(x.upper_bound),
                               float(x.objective_coefficient),
                               variable_kind_dict[x.variable_kind],
                               x.id)
                     for x in cobra_model.reactions]
    reaction_to_variable = dict(zip(cobra_model.reactions,
                                    variable_list))
    # Integrate new variables
    lp.update()
    #Set objective to quadratic program
    if quadratic_component is not None:
        if not hasattr(quadratic_component, 'todok'):
            raise Exception('quadratic component must have method todok')

        quadratic_objective = QuadExpr()
        for (index_0, index_1), the_value in quadratic_component.todok().items():
            quadratic_objective.addTerms(the_value,
                                   variable_list[index_0],
                                   variable_list[index_1])
        #Does this override the linear objective coefficients or integrate with them?
        lp.setObjective(quadratic_objective, sense=objective_sense)
    #Constraints are based on mass balance
    #Construct the lin expression lists and then add
    #TODO: Speed this up as it takes about .18 seconds
    #HERE
    for the_metabolite in cobra_model.metabolites:
        constraint_coefficients = []
        constraint_variables = []
        for the_reaction in the_metabolite._reaction:
            constraint_coefficients.append(the_reaction._metabolites[the_metabolite])
            constraint_variables.append(reaction_to_variable[the_reaction])
        #Add the metabolite to the problem
        lp.addConstr(LinExpr(constraint_coefficients, constraint_variables),
                     sense_dict[the_metabolite._constraint_sense.upper()],
                     the_metabolite._bound,
                     the_metabolite.id)



    return(lp)
###

###
def solve_problem(lp, **kwargs):
    """A performance tunable method for updating a model problem file

    """
    #Update parameter settings if provided
    if kwargs:
        [set_parameter(lp, parameter_mappings[k], v)
         for k, v in iteritems(kwargs) if k in parameter_mappings]

    try:
        print_solver_time = kwargs['print_solver_time']
        start_time = time()
    except:
        print_solver_time = False
    lp.update()
    #Different methods to try if lp_method fails
    lp.optimize()
    status = get_status(lp)
    if print_solver_time:
        print('optimize time: {:f}'.format(time() - start_time))
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
    if 'new_objective' in the_parameters and \
           the_parameters['new_objective'] not in ['update problem', None]:
       from ..flux_analysis.objective import update_objective
       update_objective(cobra_model, the_parameters['new_objective'])

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
        lp_method = 0
    the_methods = [0, 2, 1]
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
    status = solve_problem(lp, **the_parameters)
    the_solution = format_solution(lp, cobra_model)
    if status != 'optimal' and error_reporting:
        print('{:s} failed: {:s}'.format(solver_name, status))
    cobra_model.solution = the_solution
    solution = {'the_problem': lp, 'the_solution': the_solution}
    return solution
