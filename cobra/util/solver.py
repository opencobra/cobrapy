"""This module includes additional helper functions for the optlang
solver object that integrate well with the context manager, meaning that
all operations defined here are automatically reverted when used in a
`with model:` block.

The functions defined here together with the existing model functions should
allow you to implement custom flux analysis methods with ease."""

from cobra.util.context import get_context
import cobra
from functools import partial
import optlang
import re
import sympy
import six


class SolverNotFound(Exception):
    """
    A simple Exception when a solver can not be found.
    """
    pass


solvers = {match.split("_")[0]: getattr(optlang, match)
           for match in dir(optlang) if "_interface" in match}
"""
Defines all the solvers that were found in optlang.
"""


def linear_reaction_coefficients(model, reactions=None):
    """Coefficient for the reactions in a linear objective

    Parameters
    ----------
    model : cobra model
        the model object that defined the objective
    reactions : list
        an optional list for the reactions to get the coefficients for. All
        reactions if left missing.

    Returns
    -------
    dict
        A dictionary where the key is the reaction object and the value is
        the corresponding coefficient. Empty dictionary if there are no
        linear terms in the objective.
    """
    linear_coefficients = {}
    reactions = model.reactions if not reactions else reactions
    try:
        objective_expression = model.solver.objective.expression
        coefficients = objective_expression.as_coefficients_dict()
    except AttributeError:
        return linear_coefficients
    for rxn in reactions:
        forward_coefficient = coefficients.get(rxn.forward_variable, 0)
        reverse_coefficient = coefficients.get(rxn.reverse_variable, 0)
        if forward_coefficient != 0:
            if forward_coefficient == -reverse_coefficient:
                linear_coefficients[rxn] = float(forward_coefficient)
    return linear_coefficients


def set_objective(model, value, additive=False):
    """ Set the model objective

    Parameters
    ----------
    model : cobra model
       The model to set the objective for
    value : model.solver.interface.Objective,
            e.g. optlang.glpk_interface.Objective, sympy.Basic or dict

        If the model objective is linear, the value can be a new Objective
        object or a dictionary with linear coefficients where each key is a
        reaction and the element the new coefficient (float).

        If the objective is not linear and `additive` is true, only values
        of class Objective.

    additive : bool
        If true, add the terms to the current objective, otherwise start with
        an empty objective.
    """
    by_objective = isinstance(value,
                              (sympy.Basic, model.solver.interface.Objective))

    not_supported = (
        additive and (not model.objective.is_Linear and not by_objective))
    if not_supported:
        raise ValueError('can only update non-linear objectives additively '
                         'using object of class '
                         'model.solver.interface.Objective, not %s' %
                         type(value))
    reverse_value = None
    if by_objective:
        if not additive:
            if isinstance(value, sympy.Basic):
                value = model.solver.interface.Objective(value, sloppy=False)
            model.solver.objective = value
        else:
            if isinstance(value, model.solver.interface.Objective):
                value = value.expression
            model.solver.objective += value
            reverse_value = -value
    elif isinstance(value, dict):
        if not additive:
            model.solver.objective = model.solver.interface.Objective(
                sympy.S.Zero, direction='max')
        reverse_value = {}
        for reaction, coef in value.items():
            reverse_value[reaction] = reaction.objective_coefficient
            model.solver.objective.set_linear_coefficients(
                {reaction.forward_variable: coef,
                 reaction.reverse_variable: -coef})
    else:
        raise TypeError(
            '%r is not a valid objective for %r.' % (value, model.solver))
    context = get_context(model)
    if context and reverse_value:
        context(partial(set_objective, model=model, value=reverse_value,
                        additive=additive))


def interface_to_str(interface):
    """Give a string representation for an optlang interface.

    Parameters
    ----------
    interface: string
        Full name of the interface in optlang or cobra representation.
        For instance 'optlang.glpk_interface' or 'optlang-glpk'.

    Returns
    -------
    string
       The name of the interface as a string
    """
    return re.sub(r"optlang.|.interface", "", interface)


def get_solver_name(mip=False, qp=False):
    """Selects a solver for a given optimization problem in a reproducible
    manner.

    Parameters
    ----------
    mip: bool
        Does the solver require mixed integer linear programming capabilities?
    qp: bool
        Does the solver require quadratic programming capabilities?

    Returns
    -------
    string
        The name of feasible solver.

    Notes
    -----
    Raises SolverNotFound if a suitable solver is not found.
    """
    if len(solvers) == 0:
        raise SolverNotFound("no solvers installed")
    # glpk only does lp, not qp. Gurobi and cplex are better at mip
    mip_order = ["gurobi", "cplex", "mosek", "glpk"]
    lp_order = ["glpk", "cplex", "gurobi", "mosek", "scipy"]
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


def add_to_solver(model, what):
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
       `model.solver.interface.Variable` or
       `model.solver.interface.Constraint`.
    """
    context = get_context(model)

    model.solver.add(what)
    if context:
        context(partial(model.solver.remove, what))


def remove_from_solver(model, what):
    """Removes variables and constraints from a Model's solver object.

    Useful to temporarily remove variables and constraints from a Models's
    solver object.

    Parameters
    ----------
    model: a cobra model
       The model from which to remove the variables and constraints.
    what: list or tuple of optlang variables or constraints.
       The variables or constraints to remove from the model. Must be of
       class `model.solver.interface.Variable` or
       `model.solver.interface.Constraint`.
    """
    context = get_context(model)

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
       Must be a valid expression within the Model's solver object. The
       absolute value is applied automatically on the expression.
    name: string
       The name of the newly created variable.
    ub: positive float
       The upper bound for the variable.
    """
    variable = model.solver.interface.Variable(name, lb=0, ub=ub)

    # The following constraints enforce variable > expression and
    # variable > -expression
    constraints = [
        model.solver.interface.Constraint(
            expression - variable, ub=0, name="abs_pos_" + name),
        model.solver.interface.Constraint(
            expression + variable, lb=0, name="abs_neg_" + name)
    ]
    add_to_solver(model, constraints + [variable])
