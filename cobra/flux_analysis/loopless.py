# -*- coding: utf-8 -*-

"""Provides functions to remove thermodynamically infeasible loops."""

from __future__ import absolute_import

import numpy
from optlang.symbolics import Zero

from cobra.core import get_solution
from cobra.util import create_stoichiometric_matrix, nullspace


def add_loopless(model, zero_cutoff=1e-12):
    """Modify a model so all feasible flux distributions are loopless.

    In most cases you probably want to use the much faster `loopless_solution`.
    May be used in cases where you want to add complex constraints and
    objecives (for instance quadratic objectives) to the model afterwards
    or use an approximation of Gibbs free energy directions in you model.
    Adds variables and constraints to a model which will disallow flux
    distributions with loops. The used formulation is described in [1]_.
    This function *will* modify your model.

    Parameters
    ----------
    model : cobra.Model
        The model to which to add the constraints.
    zero_cutoff : positive float, optional
        Cutoff used for null space. Coefficients with an absolute value smaller
        than `zero_cutoff` are considered to be zero.

    Returns
    -------
    Nothing

    References
    ----------
    .. [1] Elimination of thermodynamically infeasible loops in steady-state
       metabolic models. Schellenberger J, Lewis NE, Palsson BO. Biophys J.
       2011 Feb 2;100(3):544-53. doi: 10.1016/j.bpj.2010.12.3707. Erratum
       in: Biophys J. 2011 Mar 2;100(5):1381.
    """
    internal = [i for i, r in enumerate(model.reactions) if not r.boundary]
    s_int = create_stoichiometric_matrix(model)[:, numpy.array(internal)]
    n_int = nullspace(s_int).T
    max_bound = max(max(abs(b) for b in r.bounds) for r in model.reactions)
    prob = model.problem

    # Add indicator variables and new constraints
    to_add = []
    for i in internal:
        rxn = model.reactions[i]
        # indicator variable a_i
        indicator = prob.Variable("indicator_" + rxn.id, type="binary")
        # -M*(1 - a_i) <= v_i <= M*a_i
        on_off_constraint = prob.Constraint(
            rxn.flux_expression - max_bound * indicator,
            lb=-max_bound, ub=0, name="on_off_" + rxn.id)
        # -(max_bound + 1) * a_i + 1 <= G_i <= -(max_bound + 1) * a_i + 1000
        delta_g = prob.Variable("delta_g_" + rxn.id)
        delta_g_range = prob.Constraint(
            delta_g + (max_bound + 1) * indicator,
            lb=1, ub=max_bound, name="delta_g_range_" + rxn.id)
        to_add.extend([indicator, on_off_constraint, delta_g, delta_g_range])

    model.add_cons_vars(to_add)

    # Add nullspace constraints for G_i
    for i, row in enumerate(n_int):
        name = "nullspace_constraint_" + str(i)
        nullspace_constraint = prob.Constraint(Zero, lb=0, ub=0, name=name)
        model.add_cons_vars([nullspace_constraint])
        coefs = {model.variables[
                 "delta_g_" + model.reactions[ridx].id]: row[i]
                 for i, ridx in enumerate(internal) if
                 abs(row[i]) > zero_cutoff}
        model.constraints[name].set_linear_coefficients(coefs)


def _add_cycle_free(model, fluxes):
    """Add constraints for CycleFreeFlux."""
    model.objective = model.solver.interface.Objective(
        Zero, direction="min", sloppy=True)
    objective_vars = []
    for rxn in model.reactions:
        flux = fluxes[rxn.id]
        if rxn.boundary:
            rxn.bounds = (flux, flux)
            continue
        if flux >= 0:
            rxn.bounds = max(0, rxn.lower_bound), max(flux, rxn.upper_bound)
            objective_vars.append(rxn.forward_variable)
        else:
            rxn.bounds = min(flux, rxn.lower_bound), min(0, rxn.upper_bound)
            objective_vars.append(rxn.reverse_variable)

    model.objective.set_linear_coefficients({v: 1.0 for v in objective_vars})


def loopless_solution(model, fluxes=None):
    """Convert an existing solution to a loopless one.

    Removes as many loops as possible (see Notes).
    Uses the method from CycleFreeFlux [1]_ and is much faster than
    `add_loopless` and should therefore be the preferred option to get loopless
    flux distributions.

    Parameters
    ----------
    model : cobra.Model
        The model to which to add the constraints.
    fluxes : dict
        A dictionary {rxn_id: flux} that assigns a flux to each reaction. If
        not None will use the provided flux values to obtain a close loopless
        solution.

    Returns
    -------
    cobra.Solution
        A solution object containing the fluxes with the least amount of
        loops possible or None if the optimization failed (usually happening
        if the flux distribution in `fluxes` is infeasible).

    Notes
    -----
    The returned flux solution has the following properties:

    - it contains the minimal number of loops possible and no loops at all if
      all flux bounds include zero
    - it has an objective value close to the original one and the same
      objective value id the objective expression can not form a cycle
      (which is usually true since it consumes metabolites)
    - it has the same exact exchange fluxes as the previous solution
    - all fluxes have the same sign (flow in the same direction) as the
      previous solution

    References
    ----------
    .. [1] CycleFreeFlux: efficient removal of thermodynamically infeasible
       loops from flux distributions. Desouki AA, Jarre F, Gelius-Dietrich
       G, Lercher MJ. Bioinformatics. 2015 Jul 1;31(13):2159-65. doi:
       10.1093/bioinformatics/btv096.
    """
    # Need to reoptimize otherwise spurious solution artifacts can cause
    # all kinds of havoc
    # TODO: check solution status
    if fluxes is None:
        sol = model.optimize(objective_sense=None)
        fluxes = sol.fluxes

    with model:
        prob = model.problem
        # Needs one fixed bound for cplex...
        loopless_obj_constraint = prob.Constraint(
            model.objective.expression,
            lb=-1e32, name="loopless_obj_constraint")
        model.add_cons_vars([loopless_obj_constraint])
        _add_cycle_free(model, fluxes)
        solution = model.optimize(objective_sense=None)
        solution.objective_value = loopless_obj_constraint.primal

    return solution


def loopless_fva_iter(model, reaction, solution=False, zero_cutoff=1e-6):
    """Plugin to get a loopless FVA solution from single FVA iteration.

    Assumes the following about `model` and `reaction`:
    1. the model objective is set to be `reaction`
    2. the model has been optimized and contains the minimum/maximum flux for
       `reaction`
    3. the model contains an auxiliary variable called "fva_old_objective"
       denoting the previous objective

    Parameters
    ----------
    model : cobra.Model
        The model to be used.
    reaction : cobra.Reaction
        The reaction currently minimized/maximized.
    solution : boolean, optional
        Whether to return the entire solution or only the minimum/maximum for
        `reaction`.
    zero_cutoff : positive float, optional
        Cutoff used for loop removal. Fluxes with an absolute value smaller
        than `zero_cutoff` are considered to be zero.

    Returns
    -------
    single float or dict
        Returns the minimized/maximized flux through `reaction` if
        all_fluxes == False (default). Otherwise returns a loopless flux
        solution containing the minimum/maximum flux for `reaction`.
    """
    current = model.objective.value
    sol = get_solution(model)
    objective_dir = model.objective.direction

    # boundary reactions can not be part of cycles
    if reaction.boundary:
        if solution:
            return sol
        else:
            return current

    with model:
        _add_cycle_free(model, sol.fluxes)
        model.slim_optimize()

        # If the previous optimum is maintained in the loopless solution it was
        # loopless and we are done
        if abs(reaction.flux - current) < zero_cutoff:
            if solution:
                return sol
            return current

        # If previous optimum was not in the loopless solution create a new
        # almost loopless solution containing only loops including the current
        # reaction. Than remove all of those loops.
        ll_sol = get_solution(model).fluxes
        reaction.bounds = (current, current)
        model.slim_optimize()
        almost_ll_sol = get_solution(model).fluxes

    with model:
        # find the reactions with loops using the current reaction and remove
        # the loops
        for rxn in model.reactions:
            rid = rxn.id
            if ((abs(ll_sol[rid]) < zero_cutoff) and
                    (abs(almost_ll_sol[rid]) > zero_cutoff)):
                rxn.bounds = max(0, rxn.lower_bound), min(0, rxn.upper_bound)

        if solution:
            best = model.optimize()
        else:
            model.slim_optimize()
            best = reaction.flux
    model.objective.direction = objective_dir
    return best
