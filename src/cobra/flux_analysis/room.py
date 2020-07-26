# -*- coding: utf-8 -*-

"""Provide regulatory on/off minimization (ROOM)."""

from __future__ import absolute_import, division

from optlang.symbolics import Zero

from cobra.flux_analysis.parsimonious import pfba


def room(model, solution=None, linear=False, delta=0.03, epsilon=1E-03):
    """
    Compute a single solution based on regulatory on/off minimization (ROOM).

    Compute a new flux distribution that minimizes the number of active
    reactions needed to accommodate a previous reference solution.
    Regulatory on/off minimization (ROOM) is generally used to assess the
    impact of knock-outs. Thus the typical usage is to provide a wildtype flux
    distribution as reference and a model in knock-out state.

    Parameters
    ----------
    model : cobra.Model
        The model state to compute a ROOM-based solution for.
    solution : cobra.Solution, optional
        A (wildtype) reference solution.
    linear : bool, optional
        Whether to use the linear ROOM formulation or not (default False).
    delta: float, optional
        The relative tolerance range (additive) (default 0.03).
    epsilon: float, optional
        The absolute tolerance range (multiplicative) (default 0.001).

    Returns
    -------
    cobra.Solution
        A flux distribution with minimal active reaction changes compared to
        the reference.

    See Also
    --------
    add_room : add ROOM constraints and objective

    """
    with model:
        add_room(model=model, solution=solution, linear=linear, delta=delta,
                 epsilon=epsilon)
        solution = model.optimize()
    return solution


def add_room(model, solution=None, linear=False, delta=0.03, epsilon=1E-03):
    r"""
    Add constraints and objective for ROOM.

    This function adds variables and constraints for applying regulatory
    on/off minimization (ROOM) to the model.

    Parameters
    ----------
    model : cobra.Model
        The model to add ROOM constraints and objective to.
    solution : cobra.Solution, optional
        A previous solution to use as a reference. If no solution is given,
        one will be computed using pFBA.
    linear : bool, optional
        Whether to use the linear ROOM formulation or not (default False).
    delta: float, optional
        The relative tolerance range which is additive in nature
        (default 0.03).
    epsilon: float, optional
        The absolute range of tolerance which is multiplicative
        (default 0.001).

    Notes
    -----
    The formulation used here is the same as stated in the original paper [1]_.
    The mathematical expression is given below:

    minimize \sum_{i=1}^m y^i
    s.t. Sv = 0
         v_min <= v <= v_max
         v_j = 0
         j ∈ A
         for 1 <= i <= m
         v_i - y_i(v_{max,i} - w_i^u) <= w_i^u      (1)
         v_i - y_i(v_{min,i} - w_i^l) <= w_i^l      (2)
         y_i ∈ {0,1}                                (3)
         w_i^u = w_i + \delta|w_i| + \epsilon
         w_i^l = w_i - \delta|w_i| - \epsilon

    So, for the linear version of the ROOM , constraint (3) is relaxed to
    0 <= y_i <= 1.

    See Also
    --------
    pfba : parsimonious FBA

    References
    ----------
    .. [1] Tomer Shlomi, Omer Berkman and Eytan Ruppin, "Regulatory on/off
     minimization of metabolic flux changes after genetic perturbations",
     PNAS 2005 102 (21) 7695-7700; doi:10.1073/pnas.0406346102

    """

    if 'room_old_objective' in model.solver.variables:
        raise ValueError('model is already adjusted for ROOM')

    # optimizes if no reference solution is provided
    if solution is None:
        solution = pfba(model)

    prob = model.problem
    variable = prob.Variable("room_old_objective", ub=solution.objective_value)
    constraint = prob.Constraint(
        model.solver.objective.expression - variable,
        ub=0.0,
        lb=0.0,
        name="room_old_objective_constraint"
    )
    model.objective = prob.Objective(Zero, direction="min", sloppy=True)
    vars_and_cons = [variable, constraint]
    obj_vars = []

    for rxn in model.reactions:
        flux = solution.fluxes[rxn.id]

        if linear:
            y = prob.Variable("y_" + rxn.id, lb=0, ub=1)
            delta = epsilon = 0.0
        else:
            y = prob.Variable("y_" + rxn.id, type="binary")

        # upper constraint
        w_u = flux + (delta * abs(flux)) + epsilon
        upper_const = prob.Constraint(
            rxn.flux_expression - y * (rxn.upper_bound - w_u),
            ub=w_u, name="room_constraint_upper_" + rxn.id)
        # lower constraint
        w_l = flux - (delta * abs(flux)) - epsilon
        lower_const = prob.Constraint(
            rxn.flux_expression - y * (rxn.lower_bound - w_l),
            lb=w_l, name="room_constraint_lower_" + rxn.id)
        vars_and_cons.extend([y, upper_const, lower_const])
        obj_vars.append(y)

    model.add_cons_vars(vars_and_cons)
    model.objective.set_linear_coefficients({v: 1.0 for v in obj_vars})
