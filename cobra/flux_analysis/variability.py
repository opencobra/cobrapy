# -*- coding: utf-8 -*-

from __future__ import absolute_import

import math
from warnings import warn
from itertools import chain

import pandas
from optlang.symbolics import Zero
from six import iteritems

from cobra.flux_analysis.loopless import loopless_fva_iter
from cobra.flux_analysis.parsimonious import add_pfba
from cobra.flux_analysis.deletion import (
    single_gene_deletion, single_reaction_deletion)
from cobra.core import get_solution
from cobra.util import solver as sutil


def flux_variability_analysis(model, reaction_list=None, loopless=False,
                              fraction_of_optimum=1.0, pfba_factor=None):
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
    if reaction_list is None:
        reaction_list = model.reactions

    fva_result = _fva_optlang(model, reaction_list, fraction_of_optimum,
                              loopless, pfba_factor)
    return pandas.DataFrame(fva_result).T


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

        m.objective = Zero  # This will trigger the reset as well
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
                           zero_cutoff=1e-9,
                           open_exchanges=False):
    """Finds reactions that cannot carry a flux with the current
    exchange reaction settings for a cobra model, using flux variability
    analysis.

    Parameters
    ----------
    model : cobra.Model
        The model to analyze
    reaction_list : list
        List of reactions to consider, use all if left missing
    zero_cutoff : float
        Flux value which is considered to effectively be zero.
    open_exchanges : bool
        If true, set bounds on exchange reactions to very high values to
        avoid that being the bottle-neck.

    Returns
    -------
    list
        List with the blocked reactions
    """
    with model:
        if open_exchanges:
            for reaction in model.exchanges:
                reaction.bounds = (min(reaction.lower_bound, -1000),
                                   max(reaction.upper_bound, 1000))
        if reaction_list is None:
            reaction_list = model.reactions
        # limit to reactions which are already 0. If the reactions already
        # carry flux in this solution, then they can not be blocked.
        model.slim_optimize()
        solution = get_solution(model, reactions=reaction_list)
        reaction_list = [rxn for rxn in reaction_list if
                         abs(solution.fluxes[rxn.id]) < zero_cutoff]
        # run fva to find reactions where both max and min are 0
        flux_span = flux_variability_analysis(
            model, fraction_of_optimum=0., reaction_list=reaction_list)
        return [rxn_id for rxn_id, min_max in flux_span.iterrows() if
                max(abs(min_max)) < zero_cutoff]


def find_essential_genes(model, threshold=None, processes=None):
    """Return a set of essential genes.

    A gene is considered essential if restricting the flux of all reactions
    that depends on it to zero causes the objective (e.g. the growth rate)
    to also be zero.

    Parameters
    ----------
    model : cobra.Model
        The model to find the essential genes for.
    threshold : float, optional
        Minimal objective flux to be considered viable. By default this is
        0.01 times the growth rate.
    processes : int, optional
        The number of parallel processes to run. Can speed up the computations
        if the number of knockouts to perform is large. If not passed,
        will be set to the number of CPUs found.

    Returns
    -------
    set
        Set of essential genes
    """
    if threshold is None:
        threshold = model.slim_optimize(error_value=None) * 1E-02
    deletions = single_gene_deletion(model, method='fba', processes=processes)
    essential = deletions.loc[deletions['growth'].isna() |
                              (deletions['growth'] < threshold), :].index
    return set(model.genes.get_by_id(g) for ids in essential for g in ids)


def find_essential_reactions(model, threshold=None, processes=None):
    """Return a set of essential reactions.

    A reaction is considered essential if restricting its flux to zero
    causes the objective (e.g. the growth rate) to also be zero.

    Parameters
    ----------
    model : cobra.Model
        The model to find the essential reactions for.
    threshold : float, optional
        Minimal objective flux to be considered viable. By default this is
        0.01 times the growth rate.
    processes : int, optional
        The number of parallel processes to run. Can speed up the computations
        if the number of knockouts to perform is large. If not passed,
        will be set to the number of CPUs found.

    Returns
    -------
    set
        Set of essential reactions
    """
    if threshold is None:
        threshold = model.slim_optimize(error_value=None) * 1E-02
    deletions = single_reaction_deletion(model, method='fba',
                                         processes=processes)
    essential = deletions.loc[deletions['growth'].isna() |
                              (deletions['growth'] < threshold), :].index
    return set(model.reactions.get_by_id(r) for ids in essential for r in ids)
