# -*- coding: utf-8 -*-

"""Provides functions to remove thermodynamically infeasible loops."""

from __future__ import absolute_import

import numpy
from six import iteritems
from sympy.core.singleton import S

from cobra.core import Metabolite, Reaction, get_solution
from cobra.util import (linear_reaction_coefficients,
                        create_stoichiometric_array)
from cobra.manipulation.modify import convert_to_irreversible
from cobra.util import nullspace


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
    s_int = create_stoichiometric_array(model)[:, numpy.array(internal)]
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
        nullspace_constraint = prob.Constraint(S.Zero, lb=0, ub=0, name=name)
        model.add_cons_vars([nullspace_constraint])
        coefs = {model.variables[
                 "delta_g_" + model.reactions[ridx].id]: row[i]
                 for i, ridx in enumerate(internal) if
                 abs(row[i]) > zero_cutoff}
        model.constraints[name].set_linear_coefficients(coefs)


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
        Note that this requires a linear objective function involving
        only the model reactions. This is the case if
        `linear_reaction_coefficients(model)` is a correct representation of
        the objective.

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
    - it has the same exact objective value as the previous solution
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
        obj_val = sol.objective_value
    else:
        coefs = linear_reaction_coefficients(model)
        obj_val = sum(coefs[rxn] * fluxes[rxn.id] for rxn in coefs)

    prob = model.problem
    with model:
        loopless_old_obj = prob.Variable("loopless_old_objective",
                                         lb=obj_val, ub=obj_val)
        loopless_obj_constraint = prob.Constraint(
            model.solver.objective.expression - loopless_old_obj,
            lb=0, ub=0, name="loopless_obj_constraint")
        model.add_cons_vars([loopless_old_obj, loopless_obj_constraint])
        model.objective = S.Zero
        for rxn in model.reactions:
            flux = fluxes[rxn.id]
            if rxn.boundary:
                rxn.bounds = (flux, flux)
                continue
            if flux >= 0:
                rxn.lower_bound = max(0, rxn.lower_bound)
                model.objective.set_linear_coefficients(
                    {rxn.forward_variable: 1, rxn.reverse_variable: -1})
            else:
                rxn.upper_bound = min(0, rxn.upper_bound)
                model.objective.set_linear_coefficients(
                    {rxn.forward_variable: -1, rxn.reverse_variable: 1})

        model.solver.objective.direction = "min"

        solution = model.optimize(objective_sense=None)

    return solution


def loopless_fva_iter(model, reaction, all_fluxes=False, zero_cutoff=1e-9):
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
    all_fluxes : boolean, optional
        Whether to return all fluxes or only the minimum/maximum for
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
    current = model.solver.objective.value

    # boundary reactions can not be part of cycles
    if reaction.boundary:
        if all_fluxes:
            return get_solution(model).fluxes
        else:
            return current

    # Reset objective to original one and get a loopless solution
    model.solver.objective.set_linear_coefficients(
        {model.variables.fva_old_objective: 1,
         reaction.forward_variable: 0, reaction.reverse_variable: 0})

    loopless = loopless_solution(model)

    # If the previous optimum is maintained in the loopless solution it was
    # loopless and we are done
    if abs(loopless[reaction.id] - current) < zero_cutoff:
        # Reset the objective to the one used in the iteration, walk around
        # the context manager for speed
        model.solver.objective.set_linear_coefficients(
            {model.variables.fva_old_objective: 0,
             reaction.forward_variable: 1, reaction.reverse_variable: -1})
        if all_fluxes:
            current = loopless
        return current

    # If previous optimum was not in the loopless solution create a new
    # almost loopless solution containing only loops including the current
    # reaction. Than remove all of those loops.
    with model:
        reaction.bounds = (current, current)
        almost_loopless = loopless_solution(model)
        # find the reactions with loops using the current reaction and remove
        # the loops
        for rxn in model.reactions:
            if ((abs(loopless[rxn.id]) < zero_cutoff) and
                    (abs(almost_loopless[rxn.id]) > zero_cutoff)):
                rxn.bounds = (0, 0)

        # Globally reset the objective to the one used in the FVA iteration
        model.solver.objective.set_linear_coefficients(
            {model.variables.fva_old_objective: 0,
             reaction.forward_variable: 1, reaction.reverse_variable: -1})

    solution = model.optimize(objective_sense=None)
    if all_fluxes:
        best = solution.fluxes
    else:
        best = reaction.flux
    return best


def construct_loopless_model(cobra_model):
    """Construct a loopless model.

    This adds MILP constraints to prevent flux from proceeding in a loop, as
    done in http://dx.doi.org/10.1016/j.bpj.2010.12.3707
    Please see the documentation for an explanation of the algorithm.

    This must be solved with an MILP capable solver.

    """
    # copy the model and make it irreversible
    model = cobra_model.copy()
    convert_to_irreversible(model)
    max_ub = max(model.reactions.list_attr("upper_bound"))
    # a dict for storing S^T
    thermo_stoic = {"thermo_var_" + metabolite.id: {}
                    for metabolite in model.metabolites}
    # Slice operator is so that we don't get newly added metabolites
    original_metabolites = model.metabolites[:]
    for reaction in model.reactions[:]:
        # Boundary reactions are not subjected to these constraints
        if len(reaction._metabolites) == 1:
            continue
        # populate the S^T dict
        bound_id = "thermo_bound_" + reaction.id
        for met, stoic in iteritems(reaction._metabolites):
            thermo_stoic["thermo_var_" + met.id][bound_id] = stoic
        # I * 1000 > v --> I * 1000 - v > 0
        reaction_ind = Reaction(reaction.id + "_indicator")
        reaction_ind.variable_kind = "integer"
        reaction_ind.upper_bound = 1
        reaction_ub = Metabolite(reaction.id + "_ind_ub")
        reaction_ub._constraint_sense = "G"
        reaction.add_metabolites({reaction_ub: -1})
        reaction_ind.add_metabolites({reaction_ub: max_ub})
        # This adds a compensating term for 0 flux reactions, so we get
        # S^T x - (1 - I) * 1001 < -1 which becomes
        # S^T x < 1000 for 0 flux reactions and
        # S^T x < -1 for reactions with nonzero flux.
        reaction_bound = Metabolite(bound_id)
        reaction_bound._constraint_sense = "L"
        reaction_bound._bound = max_ub
        reaction_ind.add_metabolites({reaction_bound: max_ub + 1})
        model.add_reaction(reaction_ind)
    for metabolite in original_metabolites:
        metabolite_var = Reaction("thermo_var_" + metabolite.id)
        metabolite_var.lower_bound = -max_ub
        model.add_reaction(metabolite_var)
        metabolite_var.add_metabolites(
            {model.metabolites.get_by_id(k): v
             for k, v in iteritems(thermo_stoic[metabolite_var.id])})
    return model
