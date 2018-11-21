# -*- coding: utf-8 -*-

"""Provide minimization of metabolic adjustment (MOMA)."""

from __future__ import absolute_import

from optlang.symbolics import Zero, add

import cobra.util.solver as sutil
from cobra.flux_analysis.parsimonious import pfba


def moma(model, solution=None, linear=True):
    """
    Compute a single solution based on (linear) MOMA.

    Compute a new flux distribution that is at a minimal distance to a
    previous reference solution. Minimization of metabolic adjustment (MOMA) is
    generally used to assess the impact
    of knock-outs. Thus the typical usage is to provide a wildtype flux
    distribution as reference and a model in knock-out state.

    Parameters
    ----------
    model : cobra.Model
        The model state to compute a MOMA-based solution for.
    solution : cobra.Solution, optional
        A (wildtype) reference solution.
    linear : bool, optional
        Whether to use the linear MOMA formulation or not (default True).

    Returns
    -------
    cobra.Solution
        A flux distribution that is at a minimal distance compared to the
        reference solution.

    See Also
    --------
    add_moma : add MOMA constraints and objective

    """
    with model:
        add_moma(model=model, solution=solution, linear=linear)
        solution = model.optimize()
    return solution


def add_moma(model, solution=None, linear=True):
    r"""Add constraints and objective representing for MOMA.

    This adds variables and constraints for the minimization of metabolic
    adjustment (MOMA) to the model.

    Parameters
    ----------
    model : cobra.Model
        The model to add MOMA constraints and objective to.
    solution : cobra.Solution, optional
        A previous solution to use as a reference. If no solution is given,
        one will be computed using pFBA.
    linear : bool, optional
        Whether to use the linear MOMA formulation or not (default True).

    Notes
    -----
    In the original MOMA [1]_ specification one looks for the flux distribution
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

    So basically we just re-center the flux space at the old solution and then
    find the flux distribution closest to the new zero (center). This is the
    same strategy as used in cameo.

    In the case of linear MOMA [2]_, we instead minimize \sum_i abs(v^t_i). The
    linear MOMA is typically significantly faster. Also quadratic MOMA tends
    to give flux distributions in which all fluxes deviate from the reference
    fluxes a little bit whereas linear MOMA tends to give flux distributions
    where the majority of fluxes are the same reference with few fluxes
    deviating a lot (typical effect of L2 norm vs L1 norm).

    The former objective function is saved in the optlang solver interface as
    ``"moma_old_objective"`` and this can be used to immediately extract the
    value of the former objective after MOMA optimization.

    See Also
    --------
    pfba : parsimonious FBA

    References
    ----------
    .. [1] Segrè, Daniel, Dennis Vitkup, and George M. Church. “Analysis of
           Optimality in Natural and Perturbed Metabolic Networks.”
           Proceedings of the National Academy of Sciences 99, no. 23
           (November 12, 2002): 15112. https://doi.org/10.1073/pnas.232349399.
    .. [2] Becker, Scott A, Adam M Feist, Monica L Mo, Gregory Hannum,
           Bernhard Ø Palsson, and Markus J Herrgard. “Quantitative
           Prediction of Cellular Metabolism with Constraint-Based Models:
           The COBRA Toolbox.” Nature Protocols 2 (March 29, 2007): 727.
    """
    if 'moma_old_objective' in model.solver.variables:
        raise ValueError('model is already adjusted for MOMA')

    # Fall back to default QP solver if current one has no QP capability
    if not linear:
        model.solver = sutil.choose_solver(model, qp=True)

    if solution is None:
        solution = pfba(model)
    prob = model.problem
    v = prob.Variable("moma_old_objective")
    c = prob.Constraint(model.solver.objective.expression - v,
                        lb=0.0, ub=0.0, name="moma_old_objective_constraint")
    to_add = [v, c]
    model.objective = prob.Objective(Zero, direction="min", sloppy=True)
    obj_vars = []
    for r in model.reactions:
        flux = solution.fluxes[r.id]
        if linear:
            components = sutil.add_absolute_expression(
                model, r.flux_expression, name="moma_dist_" + r.id,
                difference=flux, add=False)
            to_add.extend(components)
            obj_vars.append(components.variable)
        else:
            dist = prob.Variable("moma_dist_" + r.id)
            const = prob.Constraint(r.flux_expression - dist, lb=flux, ub=flux,
                                    name="moma_constraint_" + r.id)
            to_add.extend([dist, const])
            obj_vars.append(dist ** 2)
    model.add_cons_vars(to_add)
    if linear:
        model.objective.set_linear_coefficients({v: 1.0 for v in obj_vars})
    else:
        model.objective = prob.Objective(
            add(obj_vars), direction="min", sloppy=True)
