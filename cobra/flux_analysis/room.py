# -*- coding: utf-8 -*-

from __future__ import absolute_import

from optlang.symbolics import Zero


def add_room(model, solution=None, delta=0.03, epsilon=0.001, linear=False):
    '''Add constraints and objective for ROOM

    Parameters
    ----------
    model : cobra.Model
        The model to add ROOM constraints and objectve to.
    solution : cobra.Solution
        A previous solution to use as a reference.
    linear : bool
        Whether to use the linear ROOM formulation or not.

    Returns
    -------
    Nothing.

    References
    ----------
    .. [1] Tomer Shlomi, Omer Berkman and Eytan Ruppin, "Regulatory on/off
     minimization of metabolic flux changes after genetic perturbations",
     PNAS 2005 102 (21) 7695-7700; doi:10.1073/pnas.0406346102

    '''

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

        if linear is False:
            y = problem.Variable("y_" + rxn.id, type="binary")
        else:
            y = problem.Variable("y_" + rxn.id, ub=0, lb=1)
            delta = epsilon = Zero

        # upper constraint
        w_u = flux + delta * abs(flux) + epsilon
        upper_const = problem.Constraint(
            rxn.flux_expression - y * (rxn.upper_bound - w_u),
            ub=w_u,
            name="room_constraint_upper_" + rxn.id
        )
        # lower constraint
        w_l = flux - delta * abs(flux) - epsilon
        lower_const = problem.Constraint(
            rxn.flux_expression - y * (rxn.lower_bound - w_l),
            lb=w_l,
            name="room_constraint_lower_" + rxn.id
        )
        vars_and_cons.extend([y, upper_const, lower_const])
        new_objective += y

    model.add_cons_vars(vars_and_cons)
    model.objective = problem.Objective(new_objective, direction="min")
