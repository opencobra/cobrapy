# -*- coding: utf-8 -*-

from __future__ import absolute_import

import pandas
from sympy.core.singleton import S
from warnings import warn
from itertools import chain

from cobra.flux_analysis.loopless import loopless_fva_iter
from cobra.flux_analysis.parsimonious import add_pfba
from cobra.flux_analysis.single_deletion import (single_gene_deletion,
                                                 single_reaction_deletion)
from cobra.core import get_solution
from cobra.util import solver as sutil


def flux_variability_analysis(model, reaction_list=None, loopless=False,
                              fraction_of_optimum=1.0, pfba_factor=None,
                              solver=None, **solver_args):
    """Runs flux variability analysis to find the min/max flux values for each
    each reaction in `reaction_list`.

    Parameters
    ----------
    model : a cobra model
        The model for which to run the analysis. It will *not* be modified.
    reaction_list : list of cobra.Reaction or str, optional
        The reactions for which to obtain min/max fluxes. If None will use
        all reactions in the model.
    loopless : boolean, optional
        Whether to return only loopless solutions. Ignored for legacy solvers,
        also see `Notes`.
    fraction_of_optimum : float, optional
        Must be <= 1.0. Requires that the objective value is at least
        fraction * max_objective_value. A value of 0.85 for instance means that
        the objective has to be at least at 85% percent of its maximum.
    pfba_factor : float, optional
        Add additional constraint to the model that the total sum of
        absolute fluxes must not be larger than this value times the
        smallest possible sum of absolute fluxes, i.e., by setting the value
        to 1.1 then the total sum of absolute fluxes must not be more than
        10% larger than the pfba solution. Since the pfba solution is the
        one that optimally minimizes the total flux sum, the pfba_factor
        should, if set, be larger than one. Setting this value may lead to
        more realistic predictions of the effective flux bounds.
    solver : str, optional
        Name of the solver to be used. If None it will respect the solver set
        in the model (model.solver).
    **solver_args : additional arguments for legacy solver, optional
        Additional arguments passed to the legacy solver. Ignored for
        optlang solver (those can be configured using
        model.solver.configuration).

    Returns
    -------
    pandas.DataFrame
        DataFrame with reaction identifier as the index columns

        - maximum: indicating the highest possible flux
        - minimum: indicating the lowest possible flux

    Notes
    -----
    This implements the fast version as described in [1]_. Please note that
    the flux distribution containing all minimal/maximal fluxes does not have
    to be a feasible solution for the model. Fluxes are minimized/maximized
    individually and a single minimal flux might require all others to be
    suboptimal.

    Using the loopless option will lead to a significant increase in
    computation time (about a factor of 100 for large models). However, the
    algorithm used here (see [2]_) is still more than 1000x faster than the
    "naive" version using `add_loopless(model)`. Also note that if you have
    included constraints that force a loop (for instance by setting all fluxes
    in a loop to be non-zero) this loop will be included in the solution.

    References
    ----------
    .. [1] Computationally efficient flux variability analysis.
       Gudmundsson S, Thiele I.
       BMC Bioinformatics. 2010 Sep 29;11:489.
       doi: 10.1186/1471-2105-11-489, PMID: 20920235

    .. [2] CycleFreeFlux: efficient removal of thermodynamically infeasible
       loops from flux distributions.
       Desouki AA, Jarre F, Gelius-Dietrich G, Lercher MJ.
       Bioinformatics. 2015 Jul 1;31(13):2159-65.
       doi: 10.1093/bioinformatics/btv096.
    """
    legacy, solver = sutil.choose_solver(model, solver)

    if reaction_list is None:
        reaction_list = model.reactions

    if not legacy:
        fva_result = _fva_optlang(model, reaction_list, fraction_of_optimum,
                                  loopless, pfba_factor)
    else:
        if pfba_factor is not None:
            ValueError('pfba_factor only supported for optlang interfaces')
        fva_result = _fva_legacy(model, reaction_list, fraction_of_optimum,
                                 "maximize", solver, **solver_args)
    return pandas.DataFrame(fva_result).T


def _fva_legacy(cobra_model, reaction_list, fraction_of_optimum,
                objective_sense, solver, **solver_args):
    """Runs flux variability analysis to find max/min flux values

    cobra_model : :class:`~cobra.core.Model`:

    reaction_list : list of :class:`~cobra.core.Reaction`: or their id's
        The id's for which FVA should be run. If this is None, the bounds
        will be computed for all reactions in the model.

    fraction_of_optimum : fraction of optimum which must be maintained.
        The original objective reaction is constrained to be greater than
        maximal_value * fraction_of_optimum

    solver : string of solver name
        If None is given, the default solver will be used.

    """
    lp = solver.create_problem(cobra_model)
    solver.solve_problem(lp, objective_sense=objective_sense)
    solution = solver.format_solution(lp, cobra_model)
    if solution.status != "optimal":
        raise ValueError("FVA requires the solution status to be optimal, "
                         "not " + solution.status)
    # set all objective coefficients to 0
    for i, r in enumerate(cobra_model.reactions):
        if r.objective_coefficient != 0:
            f = solution.x_dict[r.id]
            new_bounds = (f * fraction_of_optimum, f)
            solver.change_variable_bounds(lp, i,
                                          min(new_bounds), max(new_bounds))
            solver.change_variable_objective(lp, i, 0.)
    return calculate_lp_variability(lp, solver, cobra_model, reaction_list,
                                    **solver_args)


def calculate_lp_variability(lp, solver, cobra_model, reaction_list,
                             **solver_args):
    """calculate max and min of selected variables in an LP"""
    fva_results = {r.id: {} for r in reaction_list}
    for what in ("minimum", "maximum"):
        sense = "minimize" if what == "minimum" else "maximize"
        for r in reaction_list:
            r_id = r.id
            i = cobra_model.reactions.index(r_id)
            solver.change_variable_objective(lp, i, 1.)
            solver.solve_problem(lp, objective_sense=sense, **solver_args)
            fva_results[r_id][what] = solver.get_objective_value(lp)
            # revert the problem to how it was before
            solver.change_variable_objective(lp, i, 0.)
    return fva_results


def _fva_optlang(model, reaction_list, fraction, loopless, pfba_factor):
    """Helper function to perform FVA with the optlang interface.

    Parameters
    ----------
    model : a cobra model
    reaction_list : list of reactions
    fraction : float, optional
        Must be <= 1.0. Requires that the objective value is at least
        fraction * max_objective_value. A value of 0.85 for instance means that
        the objective has to be at least at 85% percent of its maximum.
    loopless : boolean, optional
        Whether to return only loopless solutions.
    pfba_factor : float, optional
        Add additional constraint to the model that the total sum of
        absolute fluxes must not be larger than this value times the
        smallest possible sum of absolute fluxes, i.e., by setting the value
        to 1.1 then the total sum of absolute fluxes must not be more than
        10% larger than the pfba solution. Setting this value may lead to
        more realistic predictions of the effective flux bounds.

    Returns
    -------
    dict
        A dictionary containing the results.
    """
    prob = model.problem
    fva_results = {rxn.id: {} for rxn in reaction_list}
    with model as m:
        m.slim_optimize(error_value=None,
                        message="There is no optimal solution for the "
                                "chosen objective!")
        # Add objective as a variable to the model than set to zero
        # This also uses the fraction to create the lower bound for the
        # old objective
        fva_old_objective = prob.Variable(
            "fva_old_objective", lb=fraction * m.solver.objective.value)
        fva_old_obj_constraint = prob.Constraint(
            m.solver.objective.expression - fva_old_objective, lb=0, ub=0,
            name="fva_old_objective_constraint")
        m.add_cons_vars([fva_old_objective, fva_old_obj_constraint])

        if pfba_factor is not None:
            if pfba_factor < 1.:
                warn('pfba_factor should be larger or equal to 1', UserWarning)
            with m:
                add_pfba(m, fraction_of_optimum=0)
                ub = m.slim_optimize(error_value=None)
                flux_sum = prob.Variable("flux_sum", ub=pfba_factor * ub)
                flux_sum_constraint = prob.Constraint(
                    m.solver.objective.expression - flux_sum, lb=0, ub=0,
                    name="flux_sum_constraint")
            m.add_cons_vars([flux_sum, flux_sum_constraint])

        m.objective = S.Zero  # This will trigger the reset as well
        for what in ("minimum", "maximum"):
            sense = "min" if what == "minimum" else "max"
            for rxn in reaction_list:
                r_id = rxn.id
                rxn = m.reactions.get_by_id(r_id)
                # The previous objective assignment already triggers a reset
                # so directly update coefs here to not trigger redundant resets
                # in the history manager which can take longer than the actual
                # FVA for small models
                m.solver.objective.set_linear_coefficients(
                    {rxn.forward_variable: 1, rxn.reverse_variable: -1})
                m.solver.objective.direction = sense
                m.slim_optimize()
                sutil.check_solver_status(m.solver.status)
                if loopless:
                    value = loopless_fva_iter(m, rxn)
                else:
                    value = m.solver.objective.value
                fva_results[r_id][what] = value
                m.solver.objective.set_linear_coefficients(
                    {rxn.forward_variable: 0, rxn.reverse_variable: 0})

    return fva_results


def find_blocked_reactions(model, reaction_list=None,
                           solver=None, zero_cutoff=1e-9,
                           open_exchanges=False, **solver_args):
    """Finds reactions that cannot carry a flux with the current
    exchange reaction settings for a cobra model, using flux variability
    analysis.

    Parameters
    ----------
    model : cobra.Model
        The model to analyze
    reaction_list : list
        List of reactions to consider, use all if left missing
    solver : string
        The name of the solver to use
    zero_cutoff : float
        Flux value which is considered to effectively be zero.
    open_exchanges : bool
        If true, set bounds on exchange reactions to very high values to
        avoid that being the bottle-neck.
    **solver_args :
        Additional arguments to the solver. Ignored for optlang based solvers.

    Returns
    -------
    list
        List with the blocked reactions
    """
    legacy, solver_interface = sutil.choose_solver(model, solver)
    with model:
        if open_exchanges:
            for reaction in model.exchanges:
                reaction.bounds = (min(reaction.lower_bound, -1000),
                                   max(reaction.upper_bound, 1000))
        if reaction_list is None:
            reaction_list = model.reactions
        # limit to reactions which are already 0. If the reactions already
        # carry flux in this solution, then they can not be blocked.
        if legacy:
            solution = solver_interface.solve(model, **solver_args)
            reaction_list = [i for i in reaction_list
                             if abs(solution.x_dict[i.id]) < zero_cutoff]
        else:
            model.solver = solver_interface
            model.slim_optimize()
            solution = get_solution(model, reactions=reaction_list)
            reaction_list = [rxn for rxn in reaction_list if
                             abs(solution.fluxes[rxn.id]) < zero_cutoff]
        # run fva to find reactions where both max and min are 0
        flux_span = flux_variability_analysis(
            model, fraction_of_optimum=0., reaction_list=reaction_list,
            solver=solver, **solver_args)
        return [rxn_id for rxn_id, min_max in flux_span.iterrows() if
                max(abs(min_max)) < zero_cutoff]


def find_essential_genes(model, threshold=0.01):
    """Return a set of essential genes.

    A gene is considered essential if restricting the flux of all reactions
    that depends on it to zero causes the objective (e.g. the growth rate)
    to also be zero.

    Parameters
    ----------
    model : cobra.Model
        The model to find the essential genes for.
    threshold : float (default 0.01)
        Minimal objective flux to be considered viable.

    Returns
    -------
    set
        Set of essential genes
    """
    solution = model.optimize(raise_error=True)
    tolerance = model.solver.configuration.tolerances.feasibility
    non_zero_flux_reactions = list(
        solution[solution.fluxes.abs() > tolerance].index)
    genes_to_check = set(chain.from_iterable(
        model.reactions.get_by_id(rid).genes for rid in
        non_zero_flux_reactions))
    deletions = single_gene_deletion(model, gene_list=genes_to_check,
                                     method='fba')
    gene_ids = list(deletions[(pandas.isnull(deletions.flux)) |
                              (deletions.flux < threshold)].index)
    return set(model.genes.get_by_any(gene_ids))


def find_essential_reactions(model, threshold=0.01):
    """Return a set of essential reactions.

    A reaction is considered essential if restricting its flux to zero
    causes the objective (e.g. the growth rate) to also be zero.

    Parameters
    ----------
    model : cobra.Model
        The model to find the essential reactions for.
    threshold : float (default 0.01)
        Minimal objective flux to be considered viable.

    Returns
    -------
    set
        Set of essential reactions
    """
    solution = model.optimize(raise_error=True)
    tolerance = model.solver.configuration.tolerances.feasibility
    non_zero_flux_reactions = list(
        solution[solution.fluxes.abs() > tolerance].index)
    deletions = single_reaction_deletion(model,
                                         reaction_list=non_zero_flux_reactions,
                                         method='fba')
    reaction_ids = list(deletions[(pandas.isnull(deletions.flux)) |
                                  (deletions.flux < threshold)].index)
    return set(model.reactions.get_by_any(reaction_ids))
