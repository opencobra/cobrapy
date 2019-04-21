
"""
Find all active reactions by solving a single MILP problem
or a (small) number of LP problems

"""
import logging
from cobra.flux_analysis.loopless import fastSNP
from cobra.flux_analysis.helpers import normalize_cutoff
from optlang import Model, Variable, Constraint, Objective
from optlang.symbolics import Zero
from scipy.linalg import orth
import numpy as np

LOGGER = logging.getLogger(__name__)


def find_active_reactions(model, bigM=10000, zero_cutoff=None,
                          relax_bounds=True, solve="lp"):
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
        (default 1e-5 if relax_bounds is True, else model.tolerance).
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
        N = fastSNP(model, bigM=bigM, zero_cutoff=zero_cutoff,
                    relax_bounds=relax_bounds)
        return [model.reactions[j].id for j in range(len(model.reactions))
                if N[j, :].any()]

    max_bound = max([max(abs(r.upper_bound), abs(r.lower_bound))
                     for r in model.reactions])
    if max_bound < float("inf"):
        bigM = max(bigM, max_bound)

    eps = normalize_cutoff(model, zero_cutoff)
    try:
        # at least 100x feasibility if it is defined
        feas_tol = model.solver.configuration.tolerances.feasibility
        eps = max(eps, feas_tol * 100)
    except:
        feas_tol = eps
        LOGGER.debug("Feasibility tolerance not defined")

    if solve == "milp":
        try:
            # ensure bigM*z << eps at integrality tolerance limit
            eps = max(eps, model.solver.configuration.tolerances.integrality *
                      bigM * 10)
        except:
            LOGGER.debug("Integrality tolerance not defined")

    LOGGER.debug("parameters:\nbigM\t%.f\neps\t%.2e\nfeas_tol\t%.2e",
                 bigM, eps, feas_tol)

    with model:

        z_pos, z_neg = {}, {}
        switch_constrs = []
        prob = model.problem

        # if solving LP iteratively, fix all z+ and z- for reversible reactions
        # first to find all active irreversible reactions
        lb0 = 1 if solve == "lp" else 0
        var_type = "continuous" if solve == "lp" else "binary"

        if relax_bounds:
            cobra.flux_analysis.helpers.relax_model_bounds(model, bigM=bigM)

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

        model.add_cons_vars([z for r, z in z_pos.items()])
        model.add_cons_vars([z for r, z in z_neg.items()])
        model.add_cons_vars(switch_constrs)
        model.objective = prob.Objective(Zero, sloppy=True, direction="min")
        model.objective.set_linear_coefficients({z: 1.0 for r, z in
                                                z_pos.items()})
        model.objective.set_linear_coefficients({z: 1.0 for r, z in
                                                z_neg.items()})

        if solve == "milp":
            LOGGER.debug("Solve an MILP problem to find all active reactions")
        else:
            LOGGER.debug("Solve LP #1 to find all active irreversible" +
                         " reactions")

        sol = model.optimize()
        active_rxns = sol.fluxes[sol.fluxes.abs() >= eps *
                                 (1 - 1e-5)].index.tolist()

        LOGGER.debug("%d active reactions found", len(active_rxns))

        if solve == "lp":
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

            # to check: reversible reactions not yet found to be active
            rxns_to_check = [r for r in model.reactions if r.reversibility and
                             r.id not in active_rxns]

            # find (nearly) all forward active reversible reactions
            # un-fix their z+
            for r in rxns_to_check:
                z_pos[r].lb = 0

            sol = model.optimize()
            new_active_rxns = [r for r in rxns_to_check if r.id in
                               sol.fluxes[sol.fluxes.abs() >= eps *
                                          (1 - 1e-5)].index.tolist()]

            LOGGER.debug("Solve LP #2: min sum(z+). %d new"
                         % len(new_active_rxns) +
                         " active reversible reactions found")

            rxns_to_remove = []
            for r in rxns_to_check:
                if r in new_active_rxns:
                    # update active rxns and rxns to check
                    active_rxns.append(r.id)
                    rxns_to_remove.append(r)
                    # fix z+ and z-
                    z_pos[r].lb, z_neg[r].lb = 1, 1
                else:
                    # fix z+, un-fix z-
                    z_pos[r].lb, z_neg[r].lb = 1, 0

            for r in rxns_to_remove:
                rxns_to_check.remove(r)

            # find (nearly) all reverse active reversible reactions
            sol = model.optimize()
            new_active_rxns = [r for r in rxns_to_check if r.id in
                               sol.fluxes[sol.fluxes.abs() >= eps *
                                          (1 - 1e-5)].index.tolist()]

            LOGGER.debug("Solve LP #3: min sum(z-). %d"
                         % len(new_active_rxns) +
                         " new active reversible reactions found")

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
                    # update active rxns and rxns to check
                    active_rxns.append(r.id)
                    rxns_to_remove.append(r)

                z_pos[r].lb, z_neg[r].lb = 1, 1

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
            constr_min_flux = prob.Constraint(sum([weight_random[r] *
                                              (r.flux_expression) for r in
                                              rxns_to_check]), lb=eps)
            model.add_cons_vars(constr_min_flux)

            iter = 3

            # find any hidden forward active reversible reactions
            LOGGER.debug("Solve LPs until no forward active reversible" +
                         " reactions are found:")

            while True:
                iter += 1

                sol = model.optimize()

                if sol.status == "infeasible":
                    # no feasible solution. No any forward active reaction
                    LOGGER.debug("Solve LP #%d: infeasible." % iter +
                                 " All forward active reactions found...")

                    break

                # solution exists, update active_rxns and re-optimize
                new_active_rxns = [r for r in rxns_to_check if r.id in
                                   sol.fluxes[sol.fluxes.abs() >= feas_tol *
                                              10].index.tolist()]
                n_new = 0

                rxns_to_remove = []
                for r in new_active_rxns:
                    n_new += 1
                    # update active rxns and rxns to check
                    active_rxns.append(r.id)
                    rxns_to_remove.append(r)
                    # exclude it from the min flux constraint
                    constr_min_flux.set_linear_coefficients(
                        {r.forward_variable: 0, r.reverse_variable: 0})

                for r in rxns_to_remove:
                    rxns_to_check.remove(r)

                LOGGER.debug("Solve LP #%d:  %d new active rxns found",
                             iter, n_new)

            # find any hidden reverse active reversible reactions
            LOGGER.debug("Solve LPs until no reverse active reversible" +
                         " reactions are found:")

            # change the constraint into minimum negative flux
            constr_min_flux.lb = -eps
            constr_min_flux.ub = -eps
            constr_min_flux.lb = None

            while True:
                iter += 1

                sol = model.optimize()

                if sol.status == "infeasible":
                    # no feasible solution. No any reverse active reaction left
                    LOGGER.debug("Solve LP #%d: infeasible. Finished.", iter)
                    LOGGER.debug("%d active reactions in total.",
                                 len(active_rxns))

                    break

                # solution exists, update active_rxns and re-optimize
                new_active_rxns = [r for r in rxns_to_check if r.id in
                                   sol.fluxes[sol.fluxes.abs() >= feas_tol *
                                              1e1].index.tolist()]
                n_new = 0

                rxns_to_remove = []
                for r in new_active_rxns:
                    n_new += 1
                    # update active rxns and rxns to check
                    active_rxns.append(r.id)
                    rxns_to_remove.append(r)
                    # exclude it from the min flux constraint
                    constr_min_flux.set_linear_coefficients(
                        {r.forward_variable: 0, r.reverse_variable: 0})

                for r in rxns_to_remove:
                    rxns_to_check.remove(r)

                LOGGER.debug("Solve LP #%d: %d new active rxns found",
                             iter, n_new)

    return active_rxns


def find_reactions_in_cycles(model, bigM=10000, zero_cutoff=None,
                             relax_bounds=True, solve="lp"):
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

    Returns
    -------
    List of reactions that are in any internal cycles
    """

    with model:
        for r in model.reactions:
            if len(r.metabolites) <= 1:
                r.upper_bound, r.lower_bound = 0, 0

        return find_active_reactions(model, bigM=bigM, zero_cutoff=zero_cutoff,
                                     relax_bounds=relax_bounds, solve=solve)
