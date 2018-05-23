# -*- coding: utf-8 -*-

"""Contains functions to run Regulatory On/Off Minimization (ROOM)."""

from __future__ import absolute_import

from optlang.symbolics import Zero, NegativeOne, Add, Mul


def add_room(model, solution=None, linear=False, delta=0.03, epsilon=1E-03):
    """
    Add constraints and objective for ROOM.

    This function adds variables and constraints for applying regulatory
    on/off minimization (ROOM) to the model.

    Parameters
    ----------
    model : cobra.Model
        The model to add ROOM constraints and objectve to.
    solution : cobra.Solution, optional
        A previous solution to use as a reference.
    linear : bool, optional
        Whether to use the linear ROOM formulation or not (default False).
    delta: float, optional
        The relative tolerance range (default 0.03).
    epsilon: float, optional
        The absolute tolerance range (default 0.001).

    Notes
    -----
    The formulation used here is the same as stated in the original paper.
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
        solution = model.optimize()

    problem = model.problem
    variable = problem.Variable("room_old_objective",
                                ub=solution.objective_value)
    constraint = problem.Constraint(
        model.solver.objective.expression - variable,
        ub=0.0,
        lb=0.0,
        name="room_old_objective_constraint"
    )
    new_objective = Zero
    vars_and_cons = [variable, constraint]

    for rxn in model.reactions:
        flux = solution.fluxes[rxn.id]

        if linear:
            y = problem.Variable("y_" + rxn.id, lb=0, ub=1)
            delta = epsilon = Zero
        else:
            y = problem.Variable("y_" + rxn.id, type="binary")

        # upper constraint
        w_u = flux + (delta * abs(flux)) + epsilon
        upper_const = problem.Constraint(
            Add(rxn.flux_expression,
                Mul(NegativeOne,
                    y,
                    Add(rxn.upper_bound,
                        Mul(NegativeOne,
                            w_u)))),
            ub=w_u,
            name="room_constraint_upper_" + rxn.id
        )
        # lower constraint
        w_l = flux - (delta * abs(flux)) - epsilon
        lower_const = problem.Constraint(
            Add(rxn.flux_expression,
                Mul(NegativeOne,
                    y,
                    Add(rxn.lower_bound,
                        Mul(NegativeOne,
                            w_l)))),
            lb=w_l,
            name="room_constraint_lower_" + rxn.id,
        )
        vars_and_cons.extend([y, upper_const, lower_const])
        new_objective += y

    model.add_cons_vars(vars_and_cons)
    model.objective = problem.Objective(new_objective, direction="min")
