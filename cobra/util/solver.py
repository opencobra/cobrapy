"""This module includes additional helper functions for the optlang
solver object that integrate well with the context manager, meaning that
all operations defined here are automatically reverted when used in a
`with model:` block.

The functions defined here together with the existing model functions should
allow you to implement custom flux analysis methods with ease."""

from cobra.util.context import get_context
from functools import partial
import optlang
import re


class SolverNotFound(Exception):
    """
    A simple Exception when a solver can not be found.
    """
    None


solvers = {match.split("_")[0]: getattr(optlang, match)
           for match in dir(optlang) if "_interface" in match}
"""
Defines all the solvers that were found in optlang.
"""


def interface_to_str(interface):
    """Give a string representation for an optlang interface.

    Parameters
    ----------
    interface: string
    Full name of the interface in optlang or cobra representation.
    For instance 'optlang.glpk_interface' or 'optlang-glpk'.

    Returns
    -------
    name: string
    The name of the interface as a string
    """
    return re.sub(r"optlang.|.interface", "", interface)


def get_solver_name(mip=False, qp=False):
    """Selects a solver for a given optimization problem in a reproducible
    manner.

    Parameters
    ----------
    mip: string
    Does the solver require mixed integer linear programming capabilities?
    qp: string
    Does the solver require quadratic programming capabilities?

    Returns
    -------
    name: string
    The name of feasible solver.

    Notes
    -----
    Raises SolverNotFound if a suitable solver is not found.
    """
    if len(solvers) == 0:
        raise SolverNotFound("no solvers installed")
    # glpk only does lp, not qp. Gurobi and cplex are better at mip
    mip_order = ["gurobi", "cplex", "mosek", "glpk"]
    lp_order = ["glpk", "cplex",  "gurobi", "mosek", "scipy"]
    qp_order = ["gurobi", "cplex", "mosek"]

    if mip is False and qp is False:
        for solver_name in lp_order:
            if solver_name in solvers:
                return solver_name
        # none of them are in the list order - so return the first one
        return list(solvers)[0]
    elif qp:  # mip does not yet matter for this determination
        for solver_name in qp_order:
            if solver_name in solvers:
                return solver_name
        raise SolverNotFound("no qp-capable solver found")
    else:
        for solver_name in mip_order:
            if solver_name in solvers:
                return solver_name
    raise SolverNotFound("no mip-capable solver found")


def add_to_solver(model, what=[]):
    """Adds variables and constraints to a Model's solver object.

    Useful for variables and constraints that can not be expressed with
    reactions and lower/upper bounds. Will integrate with the Model's context
    manager in order to revert changes upon leaving the context.

    Parameters
    ----------
    model: a cobra model
    The model to which to add the variables and constraints.
    what: list or tuple of optlang variables or constraints.
    The variables or constraints to add to the model. Must be of class
    `model.solver.interface.Variable` or `model.solver.interface.Constraint`.
    """
    context = get_context(model)

    if len(what) > 0:
        model.solver.add(what)
        if context:
            context(partial(model.solver.remove, what))


def remove_from_solver(model, what=[]):
    """Removes variables and constraints from a Model's solver object.

    Useful to temporarily remove variables and constraints from a Models's
    solver object.

    Parameters
    ----------
    model: a cobra model
    The model from which to remove the variables and constraints.
    what: list or tuple of optlang variables or constraints.
    The variables or constraints to remove from the model. Must be of class
    `model.solver.interface.Variable` or `model.solver.interface.Constraint`.
    """
    context = get_context(model)

    if len(what) > 0:
        model.solver.remove(what)
        if context:
            context(partial(model.solver.add, what))


def add_absolute_expression(model, expression, name="abs_var", ub=None):
    """Adds the absolute value of an expression to the model and defines
    a variable for the absolute value that can be used in other objectives or
    constraints.

    Parameters
    ----------
    model: a cobra model
    The model to which to add the absolute expression.
    expression: A sympy expression
    Must be a valid expression within the Model's solver object. The absolute
    value is applied automatically on the expression.
    name: string
    The name of the newly created variable.
    ub: positive float
    The upper bound for the variable.
    """

    context = get_context(model)
    variable = model.solver.interface.Variable(name, lb=0, ub=ub)

    # The following constraints enforce variable > expression and
    # variable > -expression
    constraints = [
        # positive value constraint
        model.solver.interface.Constraint(expression - variable, ub=0,
                                          name="abs_pos_" + name),
        # negative value constraint
        model.solver.interface.Constraint(expression + variable, lb=0,
                                          name="abs_neg_" + name)
        ]
    model.solver.add(constraints + [variable])

    if context:
        context(partial(model.solver.remove, constraints + [variable]))
