"""Additional helper functions for the optlang solvers.

All functions integrate well with the context manager, meaning that
all operations defined here are automatically reverted when used in a
`with model:` block.

The functions defined here together with the existing model functions
should allow you to implement custom flux analysis methods with ease.

"""

import logging
import re
from functools import partial
from types import ModuleType
from typing import TYPE_CHECKING, Dict, List, NamedTuple, Optional, Tuple, Union
from warnings import warn

import optlang
import pandas as pd
from optlang.interface import (
    FEASIBLE,
    INFEASIBLE,
    ITERATION_LIMIT,
    NUMERIC,
    OPTIMAL,
    SUBOPTIMAL,
    TIME_LIMIT,
)
from optlang.symbolics import Basic, Zero

from cobra.exceptions import (
    OPTLANG_TO_EXCEPTIONS_DICT,
    OptimizationError,
    SolverNotFound,
)
from cobra.util.context import get_context


# Used to avoid cyclic reference and enable third-party static type checkers to work
if TYPE_CHECKING:
    from cobra import Model, Reaction


CONS_VARS = Union[optlang.interface.Constraint, optlang.interface.Variable]

logger = logging.getLogger(__name__)

# Define all the solvers that are found in optlang.
solvers = {
    match.split("_interface")[0]: getattr(optlang, match)
    for match in dir(optlang)
    if "_interface" in match
}

# Defines all the QP solvers implemented in optlang.
qp_solvers = ["cplex", "gurobi", "osqp"]

# optlang solution statuses which still allow retrieving primal values
has_primals = [NUMERIC, FEASIBLE, INFEASIBLE, SUBOPTIMAL, ITERATION_LIMIT, TIME_LIMIT]


class Components(NamedTuple):
    """Define an object for adding absolute expressions."""

    variable: optlang.interface.Variable
    upper_constraint: optlang.interface.Constraint
    lower_constraint: optlang.interface.Constraint


def linear_reaction_coefficients(
    model: "Model", reactions: Optional[List["Reaction"]] = None
) -> Dict["Reaction", float]:
    """Retrieve coefficient for the reactions in a linear objective.

    Parameters
    ----------
    model : cobra.Model
        The cobra model defining the linear objective.
    reactions : list of cobra.Reaction, optional
        An optional list of the reactions to get the coefficients for.
        By default, all reactions are considered (default None).

    Returns
    -------
    dict
        A dictionary where the keys are the reaction objects and the values
        are the corresponding coefficient. Empty dictionary if there are no
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


def _valid_atoms(model: "Model", expression: optlang.symbolics.Basic) -> bool:
    """Check whether a sympy expression references the correct variables.

    Parameters
    ----------
    model : cobra.Model
        The model in which to check for variables.
    expression : sympy.Basic
        A sympy expression.

    Returns
    -------
    bool
        True if all referenced variables are contained in model, False
        otherwise.

    """
    atoms = expression.atoms(optlang.interface.Variable)
    return all(a.problem is model.solver for a in atoms)


def set_objective(
    model: "Model",
    value: Union[
        optlang.interface.Objective,
        optlang.symbolics.Basic,
        Dict["Reaction", float],
    ],
    additive: bool = False,
) -> None:
    """Set the model objective.

    Parameters
    ----------
    model : cobra.Model
       The model to set the objective for.
    value : optlang.interface.Objective, optlang.symbolics.Basic, dict
        If the model objective is linear, then the value can be a new
        optlang.interface.Objective or a dictionary with linear
        coefficients where each key is a reaction and the corresponding
        value is the new coefficient (float).
        If the objective is non-linear and `additive` is True, then only
        values of class optlang.interface.Objective, are accepted.
    additive : bool
        If True, add the terms to the current objective, otherwise start with
        an empty objective.

    Raises
    ------
    ValueError
        If model objective is non-linear and the `value` is a dict.
    TypeError
        If the type of `value` is not one of the accepted ones.

    """
    interface = model.problem
    reverse_value = model.solver.objective.expression
    reverse_value = interface.Objective(
        reverse_value, direction=model.solver.objective.direction, sloppy=True
    )

    if isinstance(value, dict):
        if not model.objective.is_Linear:
            raise ValueError(
                "You can only update non-linear objectives additively using object of "
                f"class optlang.interface.Objective, not of {type(value)}"
            )

        if not additive:
            model.solver.objective = interface.Objective(
                Zero, direction=model.solver.objective.direction
            )
        for reaction, coef in value.items():
            model.solver.objective.set_linear_coefficients(
                {reaction.forward_variable: coef, reaction.reverse_variable: -coef}
            )

    elif isinstance(value, (Basic, optlang.interface.Objective)):
        if isinstance(value, Basic):
            value = interface.Objective(
                value, direction=model.solver.objective.direction, sloppy=False
            )
        # Check whether expression only uses variables from current model;
        # clone the objective if not, faster than cloning without checking
        if not _valid_atoms(model, value.expression):
            value = interface.Objective.clone(value, model=model.solver)

        if not additive:
            model.solver.objective = value
        else:
            model.solver.objective += value.expression
    else:
        raise TypeError(f"{value} is not a valid objective for {model.solver}.")

    context = get_context(model)
    if context:

        def reset():
            model.solver.objective = reverse_value
            model.solver.objective.direction = reverse_value.direction

        context(reset)


def interface_to_str(interface: Union[str, ModuleType]) -> str:
    """Give a string representation for an optlang interface.

    Parameters
    ----------
    interface : str, ModuleType
        Full name of the interface in optlang or cobra representation.
        For instance, 'optlang.glpk_interface' or 'optlang-glpk'.

    Returns
    -------
    str
       The name of the interface as a string.
    """
    if isinstance(interface, ModuleType):
        interface = interface.__name__
    return re.sub(r"optlang.|.interface", "", interface)


def get_solver_name(mip: bool = False, qp: bool = False) -> str:
    """Select a solver for a given optimization problem.

    Parameters
    ----------
    mip : bool
        True if the solver requires mixed integer linear programming capabilities.
    qp : bool
        True if the solver requires quadratic programming capabilities.

    Returns
    -------
    str
        The name of the feasible solver.

    Raises
    ------
    SolverNotFound
        If no suitable solver could be found.

    """
    if len(solvers) == 0:
        raise SolverNotFound("No solvers found.")
    # Those lists need to be updated as optlang implements more solvers
    mip_order = ["gurobi", "cplex", "glpk"]
    lp_order = ["glpk", "cplex", "gurobi"]
    qp_order = ["gurobi", "cplex", "osqp"]

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
        raise SolverNotFound("No QP-capable solver found.")
    else:
        for solver_name in mip_order:
            if solver_name in solvers:
                return solver_name
    raise SolverNotFound("No MIP-capable solver found.")


def choose_solver(
    model: "Model", solver: Optional[str] = None, qp: bool = False
) -> ModuleType:
    """Choose a solver given a solver name and model.

    This will choose a solver compatible with the model and required
    capabilities. Also respects model.solver where it can.

    Parameters
    ----------
    model : cobra.Model
        The model for which to choose the solver.
    solver : str, optional
        The name of the solver to be used (default None).
    qp : boolean, optional
        True if the solver needs quadratic programming capabilities
        (default False).

    Returns
    -------
    optlang.interface
        Valid solver for the problem.

    Raises
    ------
    SolverNotFound
        If no suitable solver could be found.

    """
    if solver is None:
        solver = model.problem
    else:
        model.solver = solver

    # Check for QP, raise error if no QP solver found
    if qp and interface_to_str(solver) not in qp_solvers:
        solver = solvers[get_solver_name(qp=True)]

    return solver


def check_solver(obj):
    """Check whether the chosen solver is valid.

    Check whether chosen solver is valid and also warn when using
    a specialized solver. Will return the optlang interface for the
    requested solver.

    Parameters
    ----------
    obj : str or optlang.interface or optlang.interface.Model
        The chosen solver.

    Raises
    ------
    SolverNotFound
        If the solver is not valid.
    """
    not_valid_interface = SolverNotFound(
        f"{obj} is not a valid solver interface. Pick one from {', '.join(solvers)}."
    )
    if isinstance(obj, str):
        try:
            interface = solvers[interface_to_str(obj)]
        except KeyError:
            raise not_valid_interface
    elif isinstance(obj, ModuleType) and hasattr(obj, "Model"):
        interface = obj
    elif isinstance(obj, optlang.interface.Model):
        interface = obj.interface
    else:
        raise not_valid_interface

    if interface_to_str(interface) in ["osqp", "coinor_cbc"]:
        logger.warning(
            "OSQP and CBC are specialized solvers for quadratic programming (QP) and "
            "mixed-integer programming (MIP) problems and may not perform well on "
            "general LP problems. So unless you intend to solve a QP or MIP problem, "
            "we recommend to change the solver back to a general purpose solver "
            "like `model.solver = 'glpk'` for instance."
        )

    return interface


def add_cons_vars_to_problem(
    model: "Model",
    what: Union[List[CONS_VARS], Tuple[CONS_VARS], Components],
    **kwargs,
) -> None:
    """Add variables and constraints to a model's solver object.

    Useful for variables and constraints that can not be expressed with
    reactions and lower/upper bounds. It will integrate with the model's
    context manager in order to revert changes upon leaving the context.

    Parameters
    ----------
    model : cobra.Model
       The model to which to add the variables and constraints.
    what : list or tuple of optlang.interface.Variable or
           optlang.interface.Constraint
       The variables and constraints to add to the model.
    **kwargs : keyword arguments
       Keyword arguments passed to solver's add() method.

    """
    model.solver.add(what, **kwargs)

    context = get_context(model)
    if context:
        context(partial(model.solver.remove, what))


def remove_cons_vars_from_problem(
    model: "Model",
    what: Union[List[CONS_VARS], Tuple[CONS_VARS], Components],
) -> None:
    """Remove variables and constraints from a model's solver object.

    Useful to temporarily remove variables and constraints from a model's
    solver object.

    Parameters
    ----------
    model : cobra.Model
       The model from which to remove the variables and constraints.
    what : list or tuple of optlang.interface.Variable or
           optlang.interface.Constraint
       The variables and constraints to remove from the model.

    """
    model.solver.remove(what)

    context = get_context(model)
    if context:
        context(partial(model.solver.add, what))


def add_absolute_expression(
    model: "Model",
    expression: str,
    name: str = "abs_var",
    ub: Optional[float] = None,
    difference: float = 0.0,
    add: bool = True,
) -> Components:
    """Add the absolute value of an expression to the model.

    Also defines a variable for the absolute value that can be used in
    other objectives or constraints.

    Parameters
    ----------
    model : cobra.Model
       The model to which to add the absolute expression.
    expression : str
       Must be a valid symbolic expression within the model's solver object.
       The absolute value is applied automatically on the expression.
    name : str, optional
       The name of the newly created variable (default "abs_var").
    ub : positive float, optional
       The upper bound for the variable (default None).
    difference : positive float, optional
        The difference between the expression and the variable
        (default 0.0).
    add : bool, optional
        Whether to add the variable to the model at once (default True).

    Returns
    -------
    Components
        A named tuple with variable and two constraints (upper_constraint,
        lower_constraint) describing the new variable and the constraints
        that assign the absolute value of the expression to it.

    """
    variable = model.problem.Variable(name, lb=0, ub=ub)
    # The following constraints enforce variable > expression and
    # variable > -expression
    upper_constraint = model.problem.Constraint(
        expression - variable, ub=difference, name="abs_pos_" + name
    )
    lower_constraint = model.problem.Constraint(
        expression + variable, lb=difference, name="abs_neg_" + name
    )
    to_add = Components(variable, upper_constraint, lower_constraint)
    if add:
        add_cons_vars_to_problem(model, to_add)
    return to_add


def fix_objective_as_constraint(
    model: "Model",
    fraction: float = 1.0,
    bound: Optional[float] = None,
    name: str = "fixed_objective_{}",
) -> float:
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
        The model to operate on.
    fraction : float, optional
        The fraction of the optimum the objective is allowed to reach
        (default 1.0).
    bound : float, optional
        The bound to use instead of fraction of maximum optimal value.
        If not None, `fraction` is ignored (default None).
    name : str, optional
        Name of the objective. May contain one "{}" placeholder which is
        filled with the name of the old objective
        (default "fixed_objective_{}").

    Returns
    -------
    float
        The value of the optimized objective * fraction

    """
    fix_objective_name = name.format(model.objective.name)
    if fix_objective_name in model.constraints:
        model.solver.remove(fix_objective_name)
    if bound is None:
        bound = model.slim_optimize(error_value=None) * fraction
    if model.objective.direction == "max":
        ub, lb = None, bound
    else:
        ub, lb = bound, None
    constraint = model.problem.Constraint(
        model.objective.expression, name=fix_objective_name, ub=ub, lb=lb
    )
    add_cons_vars_to_problem(model, constraint, sloppy=True)
    return bound


def check_solver_status(status: str = None, raise_error: bool = False) -> None:
    """Perform standard checks on a solver's status.

    Parameters
    ----------
    status: str, optional
        The status string obtained from the solver (default None).
    raise_error: bool, optional
        If True, raise error or display warning if False (default False).

    Returns
    -------
    None

    Warns
    -----
    UserWarning
        If `status` is not optimal and `raise_error` is set to True.

    Raises
    ------
    OptimizationError
        If `status` is None or is not optimal and `raise_error` is set to
        True.

    """
    if status == OPTIMAL:
        return None
    elif (status in has_primals) and not raise_error:
        warn(f"Solver status is '{status}'.", UserWarning)
    elif status is None:
        raise OptimizationError(
            "Model is not optimized yet or solver context has been switched."
        )
    else:
        raise OptimizationError(f"Solver status is '{status}'.")


def assert_optimal(model: "Model", message: str = "Optimization failed") -> None:
    """Assert model solver status is optimal.

    Do nothing if model solver status is optimal, otherwise throw
    appropriate exception depending on the status.

    Parameters
    ----------
    model : cobra.Model
        The model to check the solver status for.
    message : str, optional
        Message for the exception if solver status is not optimal
        (default "Optimization failed").

    Returns
    -------
    None

    Raises
    ------
    OptimizationError
       If solver status is not optimal.

    """
    status = model.solver.status
    if status != OPTIMAL:
        exception_cls = OPTLANG_TO_EXCEPTIONS_DICT.get(status, OptimizationError)
        raise exception_cls(f"{message} ({status}).")


def add_lp_feasibility(model: "Model") -> None:
    """Add a new objective and variables to ensure a feasible solution.

    The optimized objective will be zero for a feasible solution and
    otherwise represent the distance from feasibility (please see [1]_
    for more information).

    Parameters
    ----------
    model : cobra.Model
        The model whose feasibility is to be tested.

    Returns
    -------
    None

    References
    ----------
    .. [1] Gomez, Jose A., Kai Höffner, and Paul I. Barton.
    “DFBAlab: A Fast and Reliable MATLAB Code for Dynamic Flux Balance
    Analysis.” BMC Bioinformatics 15, no. 1 (December 18, 2014): 409.
    https://doi.org/10.1186/s12859-014-0409-8.

    """

    obj_vars = []
    prob = model.problem
    for met in model.metabolites:
        s_plus = prob.Variable("s_plus_" + met.id, lb=0)
        s_minus = prob.Variable("s_minus_" + met.id, lb=0)

        model.add_cons_vars([s_plus, s_minus])
        model.constraints[met.id].set_linear_coefficients({s_plus: 1.0, s_minus: -1.0})
        obj_vars.append(s_plus)
        obj_vars.append(s_minus)

    model.objective = prob.Objective(Zero, sloppy=True, direction="min")
    model.objective.set_linear_coefficients({v: 1.0 for v in obj_vars})


def add_lexicographic_constraints(
    model: "Model",
    objectives: List["Reaction"],
    objective_direction: Union[str, List[str]] = "max",
) -> pd.Series:
    """Successively optimize separate targets in a specific order.

    For each objective, optimize the model and set the optimal value as a
    constraint. Proceed in the order of the objectives given. Due to the
    specific order this is called lexicographic FBA [1]_. This procedure
    is useful for returning unique solutions for a set of important
    fluxes. Typically this is applied to exchange fluxes.

    Parameters
    ----------
    model : cobra.Model
        The model to be optimized.
    objectives : list of cobra.Reaction
        A list of reactions (or objectives) in the model for which unique
        fluxes are to be determined.
    objective_direction : str or list of str, optional
        The desired objective direction for each reaction (if a list) or
        the objective direction to use for all reactions (default "max").

    Returns
    -------
    pandas.Series
        A pandas Series containing the optimized fluxes for each of the
        given reactions in `objectives`.

    References
    ----------
    .. [1] Gomez, Jose A., Kai Höffner, and Paul I. Barton.
    “DFBAlab: A Fast and Reliable MATLAB Code for Dynamic Flux Balance
    Analysis.” BMC Bioinformatics 15, no. 1 (December 18, 2014): 409.
    https://doi.org/10.1186/s12859-014-0409-8.

    """

    if type(objective_direction) is not list:
        objective_direction = [objective_direction] * len(objectives)

    constraints = []
    for rxn_id, obj_dir in zip(objectives, objective_direction):
        model.objective = model.reactions.get_by_id(rxn_id)
        model.objective_direction = obj_dir
        constraints.append(fix_objective_as_constraint(model))

    return pd.Series(constraints, index=objectives)
