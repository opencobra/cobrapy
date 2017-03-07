# -*- coding: utf-8 -*-
# Interface to gurobipy

from __future__ import absolute_import

import platform
from multiprocessing import Process
from warnings import warn

from gurobipy import GRB, LinExpr, Model, QuadExpr
from six import iteritems, string_types

from ..core.solution import LegacySolution

try:
    # Import izip for python versions < 3.x
    from itertools import izip as zip
except ImportError:
    pass


def test_import():
    """Sometimes trying to import gurobipy can segfault. To prevent this from
    crashing everything, ensure it can be imported in a separate process."""
    try:
        import gurobipy
    except ImportError:
        pass

if platform.system() != "Windows":
    # https://github.com/opencobra/cobrapy/issues/207
    p = Process(target=test_import)
    p.start()
    p.join()
    if p.exitcode != 0:
        raise RuntimeError("importing gurobi causes a crash (exitcode %d)" %
                           p.exitcode)




try:
    from sympy import Basic, Number
except:
    class Basic:
        pass
    Number = Basic


def _float(value):
    if isinstance(value, Basic) and not isinstance(value, Number):
        return 0.
    else:
        return float(value)

solver_name = 'gurobi'
_SUPPORTS_MILP = True


# set solver-specific parameters
parameter_defaults = {'objective_sense': 'maximize',
                      'tolerance_optimality': 1e-6,
                      'tolerance_feasibility': 1e-6,
                      'tolerance_integer': 1e-9,
                      # This is primal simplex, default is -1 (automatic)
                      'lp_method': 0,
                      'verbose': False,
                      'log_file': ''}
parameter_mappings = {'log_file': 'LogFile',
                      'lp_method': 'Method',
                      'threads': 'Threads',
                      'objective_sense': 'ModelSense',
                      'output_verbosity': 'OutputFlag',
                      'verbose': 'OutputFlag',
                      'quadratic_precision': 'Quad',
                      'time_limit': 'TimeLimit',
                      'tolerance_feasibility': 'FeasibilityTol',
                      'tolerance_markowitz': 'MarkowitzTol',
                      'tolerance_optimality': 'OptimalityTol',
                      'iteration_limit': 'IterationLimit',
                      'tolerance_barrier': 'BarConvTol',
                      'tolerance_integer': 'IntFeasTol',
                      'MIP_gap_abs': 'MIPGapAbs',
                      'MIP_gap': 'MIPGap'}
# http://www.gurobi.com/documentation/5./6/reference-manual/method
METHODS = {"auto": -1, "primal": 0, "dual": 1, "barrier": 2,
           "concurrent": 3, "deterministic concurrent": 4}
variable_kind_dict = {'continuous': GRB.CONTINUOUS, 'integer': GRB.INTEGER}
sense_dict = {'E': GRB.EQUAL, 'L': GRB.LESS_EQUAL, 'G': GRB.GREATER_EQUAL}
objective_senses = {'maximize': GRB.MAXIMIZE, 'minimize': GRB.MINIMIZE}
status_dict = {GRB.OPTIMAL: 'optimal', GRB.INFEASIBLE: 'infeasible',
               GRB.UNBOUNDED: 'unbounded', GRB.TIME_LIMIT: 'time_limit'}


def get_status(lp):
    status = lp.status
    if status in status_dict:
        status = status_dict[status]
    else:
        status = 'failed'
    return status


def get_objective_value(lp):
    return lp.ObjVal


def format_solution(lp, cobra_model, **kwargs):
    status = get_status(lp)
    if status not in ('optimal', 'time_limit'):
        the_solution = LegacySolution(None, status=status)
    else:
        objective_value = lp.ObjVal
        x = [v.X for v in lp.getVars()]
        x_dict = {r.id: value for r, value in zip(cobra_model.reactions, x)}
        if lp.isMIP:
            y = y_dict = None  # MIP's don't have duals
        else:
            y = [c.Pi for c in lp.getConstrs()]
            y_dict = {m.id: value for m, value
                      in zip(cobra_model.metabolites, y)}
        the_solution = LegacySolution(objective_value, x=x, x_dict=x_dict, y=y,
                                y_dict=y_dict, status=status)
    return(the_solution)


def set_parameter(lp, parameter_name, parameter_value):
    if parameter_name == 'ModelSense' or parameter_name == "objective_sense":
        lp.setAttr('ModelSense', objective_senses[parameter_value])
    elif parameter_name == 'reuse_basis' and not parameter_value:
        lp.reset()
    else:
        parameter_name = parameter_mappings.get(parameter_name, parameter_name)
        if parameter_name == "Method" and isinstance(parameter_value,
                                                     string_types):
            parameter_value = METHODS[parameter_value]
        lp.setParam(parameter_name, parameter_value)


def change_variable_bounds(lp, index, lower_bound, upper_bound):
    variable = lp.getVarByName(str(index))
    variable.lb = lower_bound
    variable.ub = upper_bound


def change_variable_objective(lp, index, objective):
    variable = lp.getVarByName(str(index))
    variable.obj = objective


def change_coefficient(lp, met_index, rxn_index, value):
    met = lp.getConstrByName(str(met_index))
    rxn = lp.getVarByName(str(rxn_index))
    lp.chgCoeff(met, rxn, value)


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

    if 'copy_problem' in kwargs and kwargs['copy_problem']:
        lp = lp.copy()
    if 'reuse_basis' in kwargs and not kwargs['reuse_basis']:
        lp.reset()
    for the_variable, the_reaction in zip(lp.getVars(),
                                          cobra_model.reactions):
        the_variable.lb = float(the_reaction.lower_bound)
        the_variable.ub = float(the_reaction.upper_bound)
        the_variable.obj = float(the_reaction.objective_coefficient)


def create_problem(cobra_model, quadratic_component=None, **kwargs):
    """Solver-specific method for constructing a solver problem from
    a cobra.Model.  This can be tuned for performance using kwargs


    """
    lp = Model("")

    the_parameters = parameter_defaults
    if kwargs:
        the_parameters = parameter_defaults.copy()
        the_parameters.update(kwargs)

    # Set verbosity first to quiet infos on parameter changes
    if "verbose" in the_parameters:
        set_parameter(lp, "verbose", the_parameters["verbose"])
    for k, v in iteritems(the_parameters):
        set_parameter(lp, k, v)


    # Create variables
    #TODO:  Speed this up
    variable_list = [lp.addVar(_float(x.lower_bound),
                               _float(x.upper_bound),
                               float(x.objective_coefficient),
                               variable_kind_dict[x.variable_kind],
                               str(i))
                     for i, x in enumerate(cobra_model.reactions)]
    reaction_to_variable = dict(zip(cobra_model.reactions,
                                    variable_list))
    # Integrate new variables
    lp.update()

    #Constraints are based on mass balance
    #Construct the lin expression lists and then add
    #TODO: Speed this up as it takes about .18 seconds
    #HERE
    for i, the_metabolite in enumerate(cobra_model.metabolites):
        constraint_coefficients = []
        constraint_variables = []
        for the_reaction in the_metabolite._reaction:
            constraint_coefficients.append(_float(the_reaction._metabolites[the_metabolite]))
            constraint_variables.append(reaction_to_variable[the_reaction])
        #Add the metabolite to the problem
        lp.addConstr(LinExpr(constraint_coefficients, constraint_variables),
                     sense_dict[the_metabolite._constraint_sense.upper()],
                     the_metabolite._bound,
                     str(i))

    # Set objective to quadratic program
    if quadratic_component is not None:
        set_quadratic_objective(lp, quadratic_component)

    lp.update()
    return(lp)


def set_quadratic_objective(lp, quadratic_objective):
    if not hasattr(quadratic_objective, 'todok'):
        raise Exception('quadratic component must have method todok')
    variable_list = lp.getVars()
    linear_objective = lp.getObjective()
    # If there already was a quadratic expression set, this will be quadratic
    # and we need to extract the linear component
    if hasattr(linear_objective, "getLinExpr"):  # duck typing
        linear_objective = linear_objective.getLinExpr()
    gur_quadratic_objective = QuadExpr()
    for (index_0, index_1), the_value in quadratic_objective.todok().items():
        # gurobi does not multiply by 1/2 (only does v^T Q v)
        gur_quadratic_objective.addTerms(the_value * 0.5,
                                         variable_list[index_0],
                                         variable_list[index_1])
    # this adds to the existing quadratic objectives
    lp.setObjective(gur_quadratic_objective + linear_objective)

def solve_problem(lp, **kwargs):
    """A performance tunable method for updating a model problem file

    """
    #Update parameter settings if provided
    for k, v in iteritems(kwargs):
        set_parameter(lp, k, v)

    lp.update()
    lp.optimize()
    status = get_status(lp)
    return status


def solve(cobra_model, **kwargs):
    """

    """
    for i in ["new_objective", "update_problem", "the_problem"]:
        if i in kwargs:
            raise Exception("Option %s removed" % i)
    if 'error_reporting' in kwargs:
        warn("error_reporting deprecated")
        kwargs.pop('error_reporting')

    #Create a new problem
    lp = create_problem(cobra_model, **kwargs)

    ###Try to solve the problem using other methods if the first method doesn't work
    try:
        lp_method = kwargs['lp_method']
    except:
        lp_method = 0
    the_methods = [0, 2, 1]
    if lp_method in the_methods:
        the_methods.remove(lp_method)
    #Start with the user specified method
    the_methods.insert(0, lp_method)
    for the_method in the_methods:
        try:
            status = solve_problem(lp, lp_method=the_method)
        except:
            status = 'failed'
        if status == 'optimal':
            break

    return format_solution(lp, cobra_model)
