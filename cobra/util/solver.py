# -*- coding: utf-8 -*-

"""Additional helper functions for the optlang solvers.

All functions integrate well with the context manager, meaning that
all operations defined here are automatically reverted when used in a
`with model:` block.

The functions defined here together with the existing model functions should
allow you to implement custom flux analysis methods with ease.
"""

from __future__ import absolute_import

import re
from functools import partial
from types import ModuleType
from warnings import warn

import optlang
import sympy

from cobra.exceptions import OptimizationError
from cobra.util.context import get_context


class SolverNotFound(Exception):
    """A simple Exception when a solver can not be found."""

    pass


# Define all the solvers that are found in optlang.
solvers = {match.split("_")[0]: getattr(optlang, match)
           for match in dir(optlang) if "_interface" in match}

# Defines all the QP solvers implemented in optlang.
qp_solvers = ["cplex"]  # QP in gurobi not implemented yet


def linear_reaction_coefficients(model, reactions=None):
    """Coefficient for the reactions in a linear objective.

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


def _valid_atoms(model, expression):
    """Check whether a sympy expression references the correct variables.

    Parameters
    ----------
    model : cobra.Model
        The model in which to check for variables.
    expression : sympy.Basic
        A sympy expression.

    Returns
    -------
    boolean
        True if all referenced variables are contained in model, False
        otherwise.
    """
    atoms = expression.atoms(optlang.interface.Variable)
    return all(a.problem is model.solver for a in atoms)


def set_objective(model, value, additive=False):
    """Set the model objective.

    Parameters
    ----------
    model : cobra model
       The model to set the objective for
    value : model.problem.Objective,
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
    interface = model.problem
    reverse_value = model.solver.objective.expression
    reverse_value = interface.Objective(
        reverse_value, direction=model.solver.objective.direction,
        sloppy=True)

    if isinstance(value, dict):
        if not model.objective.is_Linear:
            raise ValueError('can only update non-linear objectives '
                             'additively using object of class '
                             'model.problem.Objective, not %s' %
                             type(value))

        if not additive:
            model.solver.objective = interface.Objective(
                sympy.S.Zero, direction=model.solver.objective.direction)
        for reaction, coef in value.items():
            model.solver.objective.set_linear_coefficients(
                {reaction.forward_variable: coef,
                 reaction.reverse_variable: -coef})

    elif isinstance(value, (sympy.Basic, optlang.interface.Objective)):
        if isinstance(value, sympy.Basic):
            value = interface.Objective(
                value, direction=model.solver.objective.direction,
                sloppy=False)
        # Check whether expression only uses variables from current model
        # clone the objective if not, faster than cloning without checking
        if not _valid_atoms(model, value.expression):
            value = interface.Objective.clone(value, model=model.solver)

        if not additive:
            model.solver.objective = value
        else:
            model.solver.objective += value.expression
    else:
        raise TypeError(
            '%r is not a valid objective for %r.' % (value, model.solver))

    context = get_context(model)
    if context:
        def reset():
            model.solver.objective = reverse_value
            model.solver.objective.direction = reverse_value.direction

        context(reset)


def interface_to_str(interface):
    """Give a string representation for an optlang interface.

    Parameters
    ----------
    interface : string
        Full name of the interface in optlang or cobra representation.
        For instance 'optlang.glpk_interface' or 'optlang-glpk'.

    Returns
    -------
    string
       The name of the interface as a string
    """
    if isinstance(interface, ModuleType):
        interface = interface.__name__
    return re.sub(r"optlang.|.interface", "", interface)


def get_solver_name(mip=False, qp=False):
    """Select a solver for a given optimization problem.

    Parameters
    ----------
    mip : bool
        Does the solver require mixed integer linear programming capabilities?
    qp : bool
        Does the solver require quadratic programming capabilities?

    Returns
    -------
    string
        The name of feasible solver.

    Raises
    ------
    SolverNotFound
        If no suitable solver could be found.
    """
    if len(solvers) == 0:
        raise SolverNotFound("no solvers installed")
    # Those lists need to be updated as optlang implements more solvers
    mip_order = ["gurobi", "cplex", "glpk"]
    lp_order = ["glpk", "cplex", "gurobi"]
    qp_order = ["cplex"]

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


def choose_solver(model, solver=None, qp=False):
    """Choose a solver given a solver name and model.

    This will choose a solver compatible with the model and required
    capabilities. Also respects model.solver where it can.

    Parameters
    ----------
    model : a cobra model
        The model for which to choose the solver.
    solver : str, optional
        The name of the solver to be used. Optlang solvers should be prefixed
        by "optlang-", for instance "optlang-glpk".
    qp : boolean, optional
        Whether the solver needs Quadratic Programming capabilities.

    Returns
    -------
    legacy : boolean
        Whether the returned solver is a legacy (old cobra solvers) version or
        an optlang solver (legacy = False).
    solver : a cobra or optlang solver interface
        Returns a valid solver for the problem. May be a cobra solver or an
        optlang interface.

    Raises
    ------
    SolverNotFound
        If no suitable solver could be found.
    """
    legacy = False
    if solver is None:
        solver = model.problem
    elif "optlang-" in solver:
        solver = interface_to_str(solver)
        solver = solvers[solver]
    else:
        legacy = True
        solver = legacy_solvers.solver_dict[solver]

    # Check for QP, raise error if no QP solver found
    # optlang only since old interface interprets None differently
    if qp and interface_to_str(solver) not in qp_solvers:
        solver = solvers[get_solver_name(qp=True)]

    return legacy, solver


def add_cons_vars_to_problem(model, what, **kwargs):
    """Add variables and constraints to a Model's solver object.

    Useful for variables and constraints that can not be expressed with
    reactions and lower/upper bounds. Will integrate with the Model's context
    manager in order to revert changes upon leaving the context.

    Parameters
    ----------
    model : a cobra model
       The model to which to add the variables and constraints.
    what : list or tuple of optlang variables or constraints.
       The variables or constraints to add to the model. Must be of class
       `model.problem.Variable` or
       `model.problem.Constraint`.
    **kwargs : keyword arguments
        passed to solver.add()
    """
    context = get_context(model)

    model.solver.add(what, **kwargs)
    if context:
        context(partial(model.solver.remove, what))


def remove_cons_vars_from_problem(model, what):
    """Remove variables and constraints from a Model's solver object.

    Useful to temporarily remove variables and constraints from a Models's
    solver object.

    Parameters
    ----------
    model : a cobra model
       The model from which to remove the variables and constraints.
    what : list or tuple of optlang variables or constraints.
       The variables or constraints to remove from the model. Must be of
       class `model.problem.Variable` or
       `model.problem.Constraint`.
    """
    context = get_context(model)

    model.solver.remove(what)
    if context:
        context(partial(model.solver.add, what))


def add_absolute_expression(model, expression, name="abs_var", ub=None):
    """Add the absolute value of an expression to the model.

    Also defines a variable for the absolute value that can be used in other
    objectives or constraints.

    Parameters
    ----------
    model : a cobra model
       The model to which to add the absolute expression.
    expression : A sympy expression
       Must be a valid expression within the Model's solver object. The
       absolute value is applied automatically on the expression.
    name : string
       The name of the newly created variable.
    ub : positive float
       The upper bound for the variable.
    """
    variable = model.problem.Variable(name, lb=0, ub=ub)

    # The following constraints enforce variable > expression and
    # variable > -expression
    constraints = [
        model.problem.Constraint(
            expression - variable, ub=0, name="abs_pos_" + name),
        model.problem.Constraint(
            expression + variable, lb=0, name="abs_neg_" + name)
    ]
    add_cons_vars_to_problem(model, constraints + [variable])


def fix_objective_as_constraint(model, fraction=1):
    """Fix current objective as an additional constraint.

    When adding constraints to a model, such as done in pFBA which
    minimizes total flux, these constraints can become too powerful,
    resulting in solutions that satisfy optimality but sacrifices too
    much for the original objective function. To avoid that, we can fix
    the current objective value as a constraint to ignore solutions that
    give a lower (or higher depending on the optimization direction)
    objective value than the original model.

    When done with the model as a context, the modification to the
    objective will be reverted when exiting that context.

    Parameters
    ----------
    model : cobra.Model
        The model to operate on
    fraction : float
        The fraction of the optimum the objective is allowed to reach.
    """
    fix_objective_name = 'Fixed_objective_{}'.format(model.objective.name)
    if fix_objective_name in model.constraints:
        model.solver.remove(fix_objective_name)
    solution = model.optimize()
    objective_bound = solution.objective_value * fraction
    if model.objective.direction == 'max':
        ub, lb = None, objective_bound
    else:
        ub, lb = objective_bound, None
    constraint = model.problem.Constraint(
        model.objective.expression,
        name=fix_objective_name, ub=ub, lb=lb)
    add_cons_vars_to_problem(model, constraint)


def check_solver_status(status):
    """Perform standard checks on a solver's status."""
    if status == "optimal":
        return
    elif status in ("non-optimal", "infeasible"):
        warn("solver status is '{}'".format(status), UserWarning)
    elif status is None:
        raise RuntimeError(
            "model was not optimized yet or solver context switched")
    else:
        raise OptimizationError("solver status is '{}'".format(status))


import cobra.solvers as legacy_solvers  # noqa
