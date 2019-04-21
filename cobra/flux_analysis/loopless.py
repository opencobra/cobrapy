# -*- coding: utf-8 -*-

"""Provides functions to remove thermodynamically infeasible loops."""

from __future__ import absolute_import

import logging

import numpy
from optlang.symbolics import Zero

from cobra.core import get_solution
from cobra.flux_analysis.helpers import normalize_cutoff
from cobra.util import create_stoichiometric_matrix, nullspace
from scipy.linalg import orth

LOGGER = logging.getLogger(__name__)


def add_loopless(model, zero_cutoff=None, method="fastSNP"):
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
        than `zero_cutoff` are considered to be zero (default model.tolerance).

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
    zero_cutoff = normalize_cutoff(model, zero_cutoff)

    if method == "original":
        internal = [i for i, r in enumerate(model.reactions) if not r.boundary]
        s_int = create_stoichiometric_matrix(model)[:, numpy.array(internal)]
        n_int = nullspace(s_int).T
    elif method == "fastSNP":
        with model:
            for r in model.reactions:
                if r.boundary:
                    r.lower_bound, r.upper_bound = 0, 0

            n_int = fastSNP(model).T

        internal = [i for i, r in enumerate(model.reactions)
                    if n_int[:, i].any()]
        n_int = n_int[:, numpy.array(internal)]

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


def loopless_fva_iter(model, reaction, solution=False, zero_cutoff=None):
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
        than `zero_cutoff` are considered to be zero (default model.tolerance).

    Returns
    -------
    single float or dict
        Returns the minimized/maximized flux through `reaction` if
        all_fluxes == False (default). Otherwise returns a loopless flux
        solution containing the minimum/maximum flux for `reaction`.
    """
    zero_cutoff = normalize_cutoff(model, zero_cutoff)

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


def fastSNP(model, bigM=1e4, zero_cutoff=None, eps=1e-3, N=None):
    """
    Find a minimal feasible sparse null space basis using fast sparse nullspace
    pursuit (Fast-SNP). Fast-SNP iteratively solves LP problems to find new
    feasible nullspace basis that lies outside the current nullspace until the
    entire feasible nullspace is found

    Parameters
    ----------
    model: cobra.Model
        cobra model. It will *not* be modified.
    bigM: float, optional
        a large constant for bounding the optimization problem, default 1e4.
    zero_cutoff: float, optional
        The cutoff to consider for zero flux (default model.tolerance).
    eps: float, optional
        The cutoff for ensuring the flux vector not lying in the current null
        nullspace i.e., the constraints w(I - P)v >= eps or <= -eps where P is
        the projection matrix of the current null space. Default 1e-3
    N: numpy.ndarray, optional
        Starting null space matrix. Default None, found by the algorithm

    Returns
    -------
    numpy.ndarray
        Null space matrix with rows corresponding to model.reactions

    Notes
    -----
    The algorithm is as follow:
    1.  N = empty matrix
    2.  P = A * A^{T} where A is an orthonormal basis for N
    3.  Solve the following two LP problems:
        min \sum_{j \in J}{|v_j|}
        s.t.   \sum_{j \in J}{S_ij * v_j} = 0   \forall i \in I
               LB_j <= v_j <= UB_j              \forall j \in J
               v_j <= |v_j|                     \forall j \in J
               -v_j <= |v_j|                    \forall j \in J
               w^{T} * (I - P) v >= eps or <= -eps (one constraint for one LP)
    4a. If at least one of the LPs is feasible, choose the solution flux vector
        v with min. non-zeros. N <- [N v]. Go to Step 2.
    4b. If infeasible, terminate and N is the minimal feasible null space.

    References
    ----------
    Saa, P. A., & Nielsen, L. K. (2016). Fast-SNP: a fast matrix pre-processing
    algorithm for efficient loopless flux optimization of metabolic models.
    Bioinformatics, 32(24), 3807-3814.

    """

    LOGGER.debug("Find minimal feasible sparse nullspace by Fast-SNP:")

    zero_cutoff = normalize_cutoff(model, zero_cutoff)
    with model:
        cobra.flux_analysis.helpers.relax_model_bounds(model, bigM=bigM)
        weight = numpy.mat(numpy.random.random(size=(1, len(model.reactions))))
        if N is None:
            wP = weight
        else:
            P_N = orth(N)
            wP = weight - numpy.matmul(numpy.matmul(weight, P_N),
                                       P_N.transpose())

        # w' (I - P'P) v >= eps / <= -eps
        constr_proj = model.problem.Constraint(0, lb=eps)
        model.add_cons_vars(constr_proj)

        # min sum(v_pos + v_neg)
        model.objective = model.problem.Objective(Zero, sloppy=True,
                                                  direction="min")
        model.objective.set_linear_coefficients(
            {r.forward_variable: 1.0 for r in model.reactions})
        model.objective.set_linear_coefficients(
            {r.reverse_variable: 1.0 for r in model.reactions})

        iter = 0
        while True:
            iter += 1
            # use w' (I - P'P) from the current null space as coefficients
            constr_proj.set_linear_coefficients(
                {model.reactions[i].forward_variable: wP[0, i] for i in
                 range(len(model.reactions))})
            constr_proj.set_linear_coefficients(
                {model.reactions[i].reverse_variable: -wP[0, i] for i in
                 range(len(model.reactions))})

            # find basis for using w' (I - P'P) v >= eps
            constr_proj.ub = bigM
            constr_proj.lb = eps
            constr_proj.ub = None
            sol = model.optimize()

            if sol.status == "optimal":
                x = sol.fluxes.to_numpy()
                x = x.reshape((len(x), 1))
                x[abs(x) < zero_cutoff] = 0
                x = x / abs(x[x != 0]).min()
            else:
                x = None

            # find basis for using w' (I - P'P) v <= -eps
            constr_proj.lb = -bigM
            constr_proj.ub = -eps
            constr_proj.lb = None
            sol = model.optimize()

            if sol.status == "optimal":
                y = sol.fluxes.to_numpy()
                y = y.reshape((len(y), 1))
                y[abs(y) < zero_cutoff] = 0
                y = y / abs(y[y != 0]).min()
            else:
                y = None

            # update N or quit
            if x is None and y is None:
                # no more feasible solution is found. Terminate.
                LOGGER.debug("iteration %d. No more feasible basis found.",
                             iter)
                LOGGER.debug("Finished")
                break
            elif x is None:
                N = y if N is None else numpy.concatenate((N, y), axis=1)
            elif y is None:
                N = x if N is None else numpy.concatenate((N, x), axis=1)
            else:
                # choose the sparsest solution
                if sum(x != 0) < sum(y != 0):
                    N = x if N is None else numpy.concatenate((N, x), axis=1)
                else:
                    N = y if N is None else numpy.concatenate((N, y), axis=1)

            LOGGER.debug("iteration %d. Feasible basis found.", iter)

            P_N = orth(N)
            wP = weight - numpy.matmul(numpy.matmul(weight, P_N),
                                       P_N.transpose())

    LOGGER.debug("The nullspace dimension is %d.", N.shape[1])

    return N
