# -*- coding: utf-8 -*-
##cobra.solvers.glpk_solver
#This script provides wrappers for pyglpk 0.3
from __future__ import absolute_import

from copy import deepcopy
from warnings import warn

from six import iteritems

from glpk import LPX

from cobra.core.solution import LegacySolution

try:
    # Import izip for python versions < 3.x
    from itertools import izip as zip
except ImportError:
    pass



solver_name = 'glpk'
_SUPPORTS_MILP = True

# solver specific parameters
variable_kind_dict = {'continuous': float, 'integer': int}
status_dict = {'opt': 'optimal', 'nofeas': 'infeasible', 'unbnd': 'unbounded'}
parameter_defaults = {
    'tolerance_feasibility': 1e-6,
    'tolerance_integer': 1e-9,
    'lp_method': 1
}

METHOD_TYPES = {"auto": 2, "primal": 1, "dual": 3}

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
        sol = LegacySolution(lp.obj.value, status=status)
        sol.x = [float(c.primal) for c in lp.cols]
        sol.x_dict = {c.name: c.primal for c in lp.cols}

        # return the duals as well as the primals for LPs
        if lp.kind == float:
            sol.y = [float(c.dual) for c in lp.rows]
            sol.y_dict = {c.name: c.dual for c in lp.rows}
        return sol

    return LegacySolution(None, status=status)

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
            raise ValueError("objective_sense should be "
                             "'maximize' or 'minimize'")
    else:
        # This will be made into an exception in the future
        warn("pyglpk parameters (other than objective_sense) are set "
             "during solve_problem")


def create_problem(cobra_model, **kwargs):
    """Solver-specific method for constructing a solver problem from
    a cobra.Model.  This can be tuned for performance using kwargs

    """
    metabolite_to_index = {r: i for i, r in enumerate(cobra_model.metabolites)}

    lp = LPX()        # Create empty problem instance
    lp.name = 'cobra'     # Assign symbolic name to problem
    lp.rows.add(len(cobra_model.metabolites))
    lp.cols.add(len(cobra_model.reactions))

    for r, the_metabolite in zip(lp.rows, cobra_model.metabolites):
        r.name = the_metabolite.id
        b = float(the_metabolite._bound)
        c = the_metabolite._constraint_sense
        if c == 'E':
            r.bounds = b, b     # Set metabolite to steady state levels
        elif c == 'L':
            r.bounds = None, b
        elif c == 'G':
            r.bounds = b, None
        else:
            raise ValueError("invalid constraint sense")

    objective_coefficients = []
    linear_constraints = []
    for c, the_reaction in zip(lp.cols, cobra_model.reactions):
        c.name = the_reaction.id
        c.kind = variable_kind_dict[the_reaction.variable_kind]
        c.bounds = the_reaction.lower_bound, the_reaction.upper_bound
        objective_coefficients.append(float(the_reaction.objective_coefficient))
        for metabolite, coefficient in iteritems(the_reaction._metabolites):
            metabolite_index = metabolite_to_index[metabolite]
            linear_constraints.append((metabolite_index, c.index, coefficient))

    #Add the new objective coefficients to the problem
    lp.obj[:] = objective_coefficients
    #Need to assign lp.matrix after constructing the whole list
    #linear_constraints.sort()  # if we wanted to be 100% deterministic
    lp.matrix = linear_constraints

    # make sure the objective sense is set in create_problem
    objective_sense = kwargs.get("objective_sense", "maximize")
    set_parameter(lp, "objective_sense", objective_sense)

    return lp


def change_variable_bounds(lp, index, lower_bound, upper_bound):
    lp.cols[index].bounds = (lower_bound, upper_bound)


def change_variable_objective(lp, index, objective):
    lp.obj[index] = objective


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

def change_coefficient(lp, met_index, rxn_index, value):
    col = lp.cols[rxn_index]
    mat = col.matrix
    for i, entry in enumerate(mat):
        if entry[0] == met_index:
            mat[i] = (met_index, value)
            col.matrix = mat
            return
    # need to append
    mat.append((met_index, value))
    col.matrix = mat


def solve_problem(lp, **kwargs):
    """A performance tunable method for updating a model problem file

    lp: a pyGLPK 0.3 problem

    For pyGLPK it is necessary to provide the following parameters, if they
    are not provided then the default settings will be used: tolerance_feasibility,
    tolerance_integer, lp_method, and objective_sense

    """
    parameters = parameter_defaults.copy()
    parameters.update(kwargs)
    if "quadratic_component" in parameters:
        if parameters.pop('quadratic_component') is not None:
            raise Exception('glpk cannot solve QPs')
    lp_args = {}  # only for lp
    extra_args = {}  # added to both lp and milp
    lp_args["tol_bnd"] = parameters.pop("tolerance_feasibility")
    lp_args["tol_dj"] = lp_args["tol_bnd"]
    method = parameters.pop("lp_method")
    if isinstance(method, int):
        lp_args["meth"] = method
    else:
        lp_args["meth"] = METHOD_TYPES[method]
    if "time_limit" in parameters:
        extra_args["tm_lim"] = int(parameters.pop("time_limit") * 1000)
    if "iteration_limit" in parameters:
        extra_args["it_lim"] = parameters.pop("iteration_limit")
    if "objective_sense" in parameters:
        set_parameter(lp, "objective_sense", parameters.pop("objective_sense"))
    tol_int = parameters.pop("tolerance_integer")
    if len(parameters) > 0:
        raise ValueError("Unknown parameters: " + ", ".join(parameters))
    lp_args.update(extra_args)
    # solve the problem
    lp.simplex(**lp_args)
    if lp.kind == int:
        # For MILPs, it is faster to solve LP then move to MILP
        lp.integer(tol_int=tol_int, **extra_args)
    return get_status(lp)


def solve(cobra_model, **kwargs):
    """Smart interface to optimization solver functions that will convert
    the cobra_model to a solver object, set the parameters, and try multiple
    methods to get an optimal solution before returning the solver object and
    a cobra.LegacySolution

    cobra_model: a cobra.Model

    returns a dict: {'the_problem': solver specific object, 'the_solution':
    cobra.solution for the optimization problem'}


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
        except Exception as e:
            status = 'failed'
        if status == 'optimal':
            break

    return format_solution(lp, cobra_model)
