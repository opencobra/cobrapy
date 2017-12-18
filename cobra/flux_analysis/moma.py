# -*- coding: utf-8 -*-

"""Contains functions to run minimization of metabolic adjustment (MOMA)."""

from __future__ import absolute_import

from optlang.symbolics import Zero

import cobra.util.solver as sutil


def add_moma(model, solution=None, linear=False):
    r"""Add constraints and objective representing for MOMA.

    This adds variables and constraints for the minimization of metabolic
    adjustment (MOMA) to the model.

    Parameters
    ----------
    model : cobra.Model
        The model to add MOMA constraints and objective to.
    solution : cobra.Solution
        A previous solution to use as a reference.
    linear : bool
        Whether to use the linear MOMA formulation or not.

    Returns
    -------
    Nothing.

    Notes
    -----
    In the original MOMA specification one looks for the flux distribution
    of the deletion (v^d) closest to the fluxes without the deletion (v).
    In math this means:

    minimize \sum_i (v^d_i - v_i)^2
    s.t. Sv^d = 0
         lb_i <= v^d_i <= ub_i

    Here, we use a variable transformation v^t := v^d_i - v_i. Substituting
    and using the fact that Sv = 0 gives:

    minimize \sum_i (v^t_i)^2
    s.t. Sv^d = 0
         v^t = v^d_i - v_i
         lb_i <= v^d_i <= ub_i

    So basically we just re-center the flux space at the old solution and than
    find the flux distribution closest to the new zero (center). This is the
    same strategy as used in cameo.

    In the case of linear MOMA, we instead minimize \sum_i abs(v^t_i). The
    linear MOMA is typically significantly faster. Also quadratic MOMA tends
    to give flux distributions in which all fluxes deviate from the reference
    fluxes a little bit whereas linear MOMA tends to give flux distributions
    where the majority of fluxes are the same reference which few fluxes
    deviating a lot (typical effect of L2 norm vs L1 norm).

    The former objective function is saved in the optlang solver interface as
    "moma_old_objective" and this can be used to immediately extract the value
    of the former objective after MOMA optimization.
    """
    if 'moma_old_objective' in model.solver.variables:
        raise ValueError('model is already adjusted for MOMA')

    # Fall back to default QP solver if current one has no QP capability
    if not linear:
        model.solver = sutil.choose_solver(model, qp=True)

    if solution is None:
        solution = model.optimize()
    prob = model.problem
    v = prob.Variable("moma_old_objective")
    c = prob.Constraint(model.solver.objective.expression - v,
                        lb=0.0, ub=0.0, name="moma_old_objective_constraint")
    to_add = [v, c]
    new_obj = Zero
    for r in model.reactions:
        flux = solution.fluxes[r.id]
        if linear:
            components = sutil.add_absolute_expression(
                model, r.flux_expression, name="moma_dist_" + r.id,
                difference=flux, add=False)
            to_add.extend(components)
            new_obj += components.variable
        else:
            dist = prob.Variable("moma_dist_" + r.id)
            const = prob.Constraint(r.flux_expression - dist, lb=flux, ub=flux,
                                    name="moma_constraint_" + r.id)
            to_add.extend([dist, const])
            new_obj += dist**2
    model.add_cons_vars(to_add)
    model.objective = prob.Objective(new_obj, direction='min')
