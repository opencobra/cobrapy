
"""
Find all active reactions by solving a single MILP problem
or a (small) number of LP problems

"""
import logging
from cobra.flux_analysis.loopless import fastSNP
from cobra.flux_analysis.helpers import normalize_cutoff, relax_model_bounds
from cobra.exceptions import OptimizationError
from optlang import Model, Variable, Constraint, Objective
from optlang.symbolics import Zero
from scipy.linalg import orth
import numpy as np

LOGGER = logging.getLogger(__name__)


def find_active_reactions(model, bigM=10000, zero_cutoff=None,
                          relax_bounds=True, solve="lp", max_iterations=10000):
    """
    Find all active reactions by solving a single MILP problem
    or a (small) number of LP problems

    Parameters
    ----------
    model: cobra.Model
        cobra model. It will *not* be modified.
    bigM: float, optional
        a large constant for bounding the optimization problem, default 1e4.
    zero_cutoff: float, optional
        The cutoff to consider for zero flux
        Default model.tolerance.
    relax_bounds: True or False
        Whether to relax the model bounds.
        All +ve LB set as 0, all -ve UB set as 0, all -ve LB set as -bigM,
        all +ve UB set as bigM. If False, use the original bounds in the model.
        Default True
    solve: "lp", "milp" or "fastSNP"
        - "lp":      to solve a number of LP problems
                     (usually finish by solving 5 LPs)
        - "milp":    to solve a single MILP problem
        - "fastSNP": find a minimal nullspace basis using Fast-SNP and then
                     return the reactions with nonzero entries in the nullspace
        Default "lp". It is recommanded to tighten model.tolerance (e.g. 1e-9)
        when calling this function especially for solving MILP.
    max_iterations: integer, optional
        The maximum number of iterations when solving LPs iteratively.
        Used only if solve == "lp". Default 10000

    Returns
    -------
    list
        List of reaction IDs which can carry flux.

    Notes
    -----
    The optmization problem solved is as follow:
    MILP version:
    .. math::

        min \sum_{j \in J}{z_pos_j + z_neg_j}
        s.t.
        \sum_{j \in J}{S_{ij} * v_{j}} = 0  \forall i \in I
        LB_j <= v_j <= UB_j             \forall j \in J
        v_j + \varepsilon z^{+}_j >=  \varepsilon  for j \in J with LB_j >= 0
        v_j - \varepsilon z^{-}_j <= -\varepsilon  for j \in J with UB_j <= 0
        v_j + Mz^{+}_j >=  \varepsilon  for j \in J with LB_j < 0 and UB_j > 0
        v_j - Mz^{-}_j <= -\varepsilon  for j \in J with LB_j < 0 and UB_j > 0
        v_j \in \mathbb{R}
        z^{+}_j \in \mathbb{R} for j \in J with LB_j >= 0
        z^{-}_j \in \mathbb{R} for j \in J with UB_j <= 0
        z^{+}_j, z^{-}_j \in {0,1} for j \in J with LB_j < 0 and UB_j > 0

    LP version:
    Solve a number of versions of the LP relaxation of the above problem
    as follows (cumulative changes in each step):
    1. Fix all :math:`z^{+}_j, z^{-}_j = 1` for all reversible reactions.
       Solve the LP to find all active irreversible reactions and
       some active reversible reactions.
    2. Fix all :math:`z^{+}_j, z^{-}_j = 1` for all irreversible reactions.
       Unfix :math:`z^{+}_j` for the reversible reactions not yet found active.
       Solve the LP to find some active reversible reactions
    3. Fix all math:`z^{+}_j`. Unfix :math:`z^{-}_j` for reversible reactions
       not yet found active. Solve the LP to find some active rev. reactions
    4. Add a randomly weighted min. flux constraint:
       \sum_{j \in J^{rev,not\: active}}{w_j v_j} >= \varepsilon
       Solve and update the set of reversible reactions not yet found active
       (if any) until infeasibility
    5. Change the sense and R.H.S. of the min. flux constraint in Step 4 to
       '<= -\varepsilon'. Solve and update until infeasibility

    References
    ----------
    Chan, S. H., Wang, L., Dash, S., & Maranas, C. D. (2018). Accelerating flux
    balance calculations in genome-scale metabolic models by localizing the
    application of loopless constraints. Bioinformatics, 34(24), 4248-4255.

    """

    solve = solve.lower()
    if solve not in ["lp", "milp", "fastsnp"]:
        raise ValueError("Parameter solve must be 'lp', 'milp' or 'fastSNP'.")

    elif solve == "fastsnp":
        N = fastSNP(model, bigM=bigM, zero_cutoff=zero_cutoff)
        return [model.reactions[j].id for j in range(len(model.reactions))
                if N[j, :].any()]

    with model:

        # construct the optimization problem
        z_pos, z_neg, eps = build_opt_problem(model, bigM=bigM,
                                              zero_cutoff=zero_cutoff,
                                              relax_bounds=relax_bounds,
                                              solve=solve)
        active_rxns = []
        if solve == "milp":
            LOGGER.debug("Solve an MILP problem to find all active reactions")
        else:
            LOGGER.debug("Solve LP #1 to find all active irreversible" +
                         " reactions")

        new_active_rxns = optimize_and_get_rxns(model, eps, model.reactions)
        active_rxns += [r.id for r in new_active_rxns]

        if solve == "milp":
            # finished if solving MILP
            return active_rxns

        # continue the iterative LP solution procedure
        # to check: reversible reactions not yet found to be active
        rxns_to_check = [r for r in model.reactions if r.reversibility and
                         r.id not in active_rxns]

        # find forward active reversible reactions
        setup_to_find_fwd_active_rxns(model, z_pos, z_neg, active_rxns,
                                      rxns_to_check)
        LOGGER.debug("Solve LP #2: min sum(z+) to find forward active" +
                     " reversible reactions.")
        new_active_rxns = optimize_and_get_rxns(model, eps, rxns_to_check)
        active_rxns += [r.id for r in new_active_rxns]

        # find reverse active reversible reactions
        setup_to_find_rev_active_rxns(z_pos, z_neg, new_active_rxns,
                                      rxns_to_check)
        LOGGER.debug("Solve LP #3: min sum(z-) to find reverse active" +
                     " reversible reactions.")
        new_active_rxns = optimize_and_get_rxns(model, eps, rxns_to_check)
        active_rxns += [r.id for r in new_active_rxns]

        # loop to find any hidden forward active reactions
        constr_min_flux, n_lp_solved, min_flux = loop_to_find_fwd_active_rxns(
            model, z_pos, z_neg, active_rxns, new_active_rxns, rxns_to_check,
            eps, max_iterations)

        # loop to find any hidden reverse active reactions
        loop_to_find_rev_active_rxns(model, active_rxns, rxns_to_check,
                                     constr_min_flux, n_lp_solved,
                                     eps, min_flux, max_iterations)

        return active_rxns


def find_reactions_in_cycles(model, bigM=10000, zero_cutoff=None,
                             relax_bounds=True, solve="lp",
                             max_iterations=10000):
    """
    Find all reactions that participate in any internal cycles. It is done by
    shutting down all exchange reactions run `find_active_reactions`

    Parameters
    ----------
    model: cobra.Model
        cobra model. It will *not* be modified.
    bigM: float, optional
        a large constant for bounding the optimization problem, default 1e4.
    zero_cutoff: float, optional
        The cutoff to consider for zero flux (default model.tolerance).
    relax_bounds: True or False
        Whether to relax the model bounds.
        All +ve LB set as 0, all -ve UB set as 0, all -ve LB set as -bigM,
        all +ve UB set as bigM. If False, use the original bounds in the model.
        Default True
    solve: "lp", "milp" or "fastSNP"
        - "lp":      to solve a number of LP problems
                     (usually finish by solving 5 LPs)
        - "milp":    to solve a single MILP problem
        - "fastSNP": find a minimal nullspace basis using Fast-SNP and then
                     return the reactions with nonzero entries in the nullspace
        Default "lp". For models with numerical difficulties when using "lp"
        or "milp", it is recommanded to tighten all tolerances:
        feasbility, optimality and integrality
    max_iterations: integer, optional
        max iterations for running the loops to solve LPs to find active
        reversible reactions. Used only if solve == "lp". Default 10000

    Returns
    -------
    List of reactions that are in any internal cycles
    """

    with model:
        for r in model.reactions:
            if r.boundary:
                r.upper_bound, r.lower_bound = 0, 0

        return find_active_reactions(model, bigM=bigM, zero_cutoff=zero_cutoff,
                                     relax_bounds=relax_bounds, solve=solve,
                                     max_iterations=max_iterations)


def optimize_and_get_rxns(model, eps, rxns_to_check):
    """
    Subroutine of `find_active_reactions`.
    Solve the optimization problem and get newly found active reactions.
    """

    new_active_rxns = None
    try:
        sol = model.optimize()
        if sol.status != "infeasible":
            new_active_rxns = [r for r in rxns_to_check if r.id in
                               sol.fluxes[sol.fluxes.abs() >= eps *
                                          (1 - 1e-5)].index.tolist()]
            LOGGER.debug("%d active reactions found", len(new_active_rxns))
    except OptimizationError:
        LOGGER.debug("Optimization error. Treat as infeasibility")

    return new_active_rxns


def build_opt_problem(model, bigM=10000, zero_cutoff=None, relax_bounds=True,
                      solve="lp"):
    """
    Subroutine of `find_active_reactions`.
    Build the optimization problem.
    """

    # make sure bigM is not smaller than the largest bound
    max_bound = max([max(abs(r.upper_bound), abs(r.lower_bound))
                     for r in model.reactions])
    if max_bound < float("inf"):
        bigM = max(bigM, max_bound)

    # ensure bigM*z << eps at tolerance limit
    eps = model.tolerance * bigM * 10
    if zero_cutoff is not None:
        eps = max(zero_cutoff, eps)

    LOGGER.debug("parameters:\nbigM\t%.f\neps\t%.2e\nfeas_tol\t%.2e",
                 bigM, eps, model.tolerance)

    z_pos, z_neg = {}, {}
    switch_constrs = []
    prob = model.problem

    # if solving LP iteratively, fix all z+ and z- for reversible reactions
    # first to find all active irreversible reactions
    lb0 = 1 if solve == "lp" else 0
    var_type = "continuous" if solve == "lp" else "binary"

    if relax_bounds:
        relax_model_bounds(model, bigM=bigM)

    for r in model.reactions:

        if r.upper_bound > 0 and r.lower_bound < 0:
            z_pos[r] = prob.Variable("z_pos_" + r.id, lb=lb0, ub=1,
                                     type=var_type)
            z_neg[r] = prob.Variable("z_neg_" + r.id, lb=lb0, ub=1,
                                     type=var_type)
            coeff_pos, coeff_neg = bigM, -bigM
        elif r.upper_bound > 0:
            z_pos[r] = prob.Variable("z_pos_" + r.id, lb=0, ub=1,
                                     type="continuous")
            coeff_pos = eps
        elif r.lower_bound < 0:
            z_neg[r] = prob.Variable("z_neg_" + r.id, lb=0, ub=1,
                                     type="continuous")
            coeff_neg = -eps

        if r.upper_bound > 0:
            # v +  eps * z_pos >= eps       for irreversible reactions
            # v + bigM * z_pos >= eps       for reversible reactions
            switch_constrs.append(prob.Constraint(r.flux_expression +
                                  coeff_pos * z_pos[r], lb=eps))

        if r.lower_bound < 0:
            # v -  eps * z_neg <= -eps       for irreversible reactions
            # v - bigM * z_neg <= -eps       for reversible reactions
            switch_constrs.append(prob.Constraint(r.flux_expression +
                                  coeff_neg * z_neg[r], ub=-eps))

    model.add_cons_vars([z for z in z_pos.values()])
    model.add_cons_vars([z for z in z_neg.values()])
    model.add_cons_vars(switch_constrs)
    model.objective = prob.Objective(Zero, sloppy=True, direction="min")
    model.objective.set_linear_coefficients({z: 1.0 for z in z_pos.values()})
    model.objective.set_linear_coefficients({z: 1.0 for z in z_neg.values()})

    return (z_pos, z_neg, eps)


def setup_to_find_fwd_active_rxns(model, z_pos, z_neg, active_rxns,
                                  rxns_to_check):
    """
    Subroutine of `find_active_reactions`.
    Change bounds to find (usually all) forward active reversible reactions.
    """

    for r in model.reactions:
        # fix z+ and z- for all irreversible reactions.
        # They are determined at this point
        if r.lower_bound >= 0 and r.upper_bound > 0:
            z_pos[r].lb = 1

        if r.lower_bound < 0 and r.upper_bound <= 0:
            z_neg[r].lb = 1

        # fix z+ and z- for reversible reactions found active.
        if r.reversibility and r.id in active_rxns:
            z_pos[r].lb, z_neg[r].lb = 1, 1

        # also fix the inactive irreversible reactions
        if not r.reversibility and r.id not in active_rxns:
            r.upper_bound, r.lower_bound = 0, 0

    # find (nearly) all forward active reversible reactions
    # un-fix their z+
    for r in rxns_to_check:
        z_pos[r].lb = 0


def setup_to_find_rev_active_rxns(z_pos, z_neg, new_active_rxns,
                                  rxns_to_check):
    """
    Subroutine of `find_active_reactions`.
    Change bounds to find (usually all) reverse active reversible reactions.
    """

    rxns_to_remove = []
    for r in rxns_to_check:
        if r in new_active_rxns:
            rxns_to_remove.append(r)
            # Already found active. Fix z+ and z-
            z_pos[r].lb, z_neg[r].lb = 1, 1
        else:
            # Not yet found active. Fix z+, un-fix z-
            z_pos[r].lb, z_neg[r].lb = 1, 0

    # update rxns_to_check
    for r in rxns_to_remove:
        rxns_to_check.remove(r)


def loop_to_find_fwd_active_rxns(model, z_pos, z_neg, active_rxns,
                                 new_active_rxns, rxns_to_check, eps,
                                 max_iterations):
    """
    Subroutine of `find_active_reactions`.
    Find flux distributions with >=1 positive flux until infeasibility.
    """

    # Usually all active reactions would have been found at this point.
    # But if there exists some reversible reactions such that
    # (i)  no irreversible reaction directionally coupled to them
    #     (no irr. rxn s.t. v_irr != 0 => v_rev != 0,
    #      otherwise it will be found by the 1st LP), or
    # (ii) any flux distributions with nonzero flux through these
    #      reactions must have \sum_{r \in rxns_to_check}{v_r} = 0
    #      (otherwise the previous two LPs would have identified them)
    # then it is possible that there are hidden active reactions,
    # though intuitively this should be uncommon for genome-scale
    # metabolic models.
    # Solve additional (>=2) LPs to find all hidden active reactions

    # fix all z+ and z-. Objective function is useful at this point
    rxns_to_remove = []
    for r in rxns_to_check:
        if r in new_active_rxns:
            rxns_to_remove.append(r)

        z_pos[r].lb, z_neg[r].lb = 1, 1

    # update rxns to check
    for r in rxns_to_remove:
        rxns_to_check.remove(r)

    # add a randomly weighted constraint for minimum flux
    # Remark:
    # Theoretically, there is a zero probability of getting a weight
    # vector w lying in the row space of the stoichiometric matrix
    # assuming the row space is not R^n for a network with n reactions.
    # If this happens, one cannot eliminate the possibility of false
    # negative results when the LPs below return infeasibility,
    # i.e., there is a non-zero flux vector v involving rxns_to_check
    # but <w, v> = 0. In practice, the probability of this happening is
    # very ... very tiny, though nonzero because the random numbers
    # drawn have a definite finite number of digits.
    # Otherwise, one needs to either run FVA for each of the remaining
    # reactions or solve MILP problems, e.g., the MILP problem above
    # or `minSpan` in Bordbar, A. et al. (2014) Minimal metabolic
    # pathway structure is consistent with associated biomolecular
    # interactions. Molecular systems biology, 10(7), 737.)

    weight_random = {r: np.round(np.random.random(), 6) for r in
                     rxns_to_check}
    min_flux = sum(w for w in weight_random.values()) * model.tolerance * 10
    min_flux = max(eps, min_flux)
    LOGGER.debug("min_flux: %.4e" % min_flux)
    constr_min_flux = model.problem.Constraint(
                      sum(w * r.flux_expression for r,w in
                          weight_random.items()), lb=min_flux)
    model.add_cons_vars(constr_min_flux)

    n_lp_solved = 3
    feas_tol = model.tolerance

    # find any hidden forward active reversible reactions
    LOGGER.debug("Solve LPs until all forward active reversible" +
                 " reactions are found:")

    for i in range(max_iterations):
        n_lp_solved += 1
        LOGGER.debug("Solve LP #%d:" % n_lp_solved)
        new_active_rxns = optimize_and_get_rxns(model, eps, rxns_to_check)

        if new_active_rxns is None:
            # no feasible solution. No any forward active reaction
            LOGGER.debug("All forward active reactions found...")

            break

        # solution exists, update active_rxns and re-optimize
        for r in new_active_rxns:
            # update active rxns and rxns to check
            active_rxns.append(r.id)
            rxns_to_check.remove(r)
            # exclude it from the min flux constraint
            constr_min_flux.set_linear_coefficients(
                {r.forward_variable: 0, r.reverse_variable: 0})

    return (constr_min_flux, n_lp_solved, min_flux)


def loop_to_find_rev_active_rxns(model, active_rxns, rxns_to_check,
                                 constr_min_flux, n_lp_solved,
                                 eps, min_flux, max_iterations):
    """
    Subroutine of `find_active_reactions`.
    Find flux distributions with >=1 negative flux until infeasibility.
    """

    # find any hidden reverse active reversible reactions
    LOGGER.debug("Solve LPs until all reverse active reversible" +
                 " reactions are found:")

    # change the constraint into minimum negative flux
    constr_min_flux.lb = -min_flux
    constr_min_flux.ub = -min_flux
    constr_min_flux.lb = None
    feas_tol = model.tolerance

    for i in range(max_iterations):
        n_lp_solved += 1
        LOGGER.debug("Solve LP #%d:" % n_lp_solved)
        new_active_rxns = optimize_and_get_rxns(model, eps, rxns_to_check)

        if new_active_rxns is None:
            # no feasible solution. No any reverse active reaction left
            LOGGER.debug("All forward active reactions found...")
            LOGGER.debug("%d active reactions in total.",
                         len(active_rxns))
            break

        # solution exists, update active_rxns and re-optimize
        for r in new_active_rxns:
            # update active rxns and rxns to check
            active_rxns.append(r.id)
            rxns_to_check.remove(r)
            # exclude it from the min flux constraint
            constr_min_flux.set_linear_coefficients(
                {r.forward_variable: 0, r.reverse_variable: 0})
