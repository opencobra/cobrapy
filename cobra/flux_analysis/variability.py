# -*- coding: utf-8 -*-

from __future__ import absolute_import

import logging
import multiprocessing
from builtins import map
from warnings import warn

from numpy import zeros
from optlang.symbolics import Zero
from pandas import DataFrame

from cobra.core import Configuration, get_solution
from cobra.flux_analysis.deletion import (
    single_gene_deletion, single_reaction_deletion)
from cobra.flux_analysis.helpers import normalize_cutoff
from cobra.flux_analysis.loopless import loopless_fva_iter
from cobra.flux_analysis.parsimonious import add_pfba
from cobra.util import solver as sutil


LOGGER = logging.getLogger(__name__)
CONFIGURATION = Configuration()


def _init_worker(model, loopless, sense):
    """Initialize a global model object for multiprocessing."""
    global _model
    global _loopless
    _model = model
    _model.solver.objective.direction = sense
    _loopless = loopless


def _fva_step(reaction_id):
    global _model
    global _loopless
    rxn = _model.reactions.get_by_id(reaction_id)
    # The previous objective assignment already triggers a reset
    # so directly update coefs here to not trigger redundant resets
    # in the history manager which can take longer than the actual
    # FVA for small models
    _model.solver.objective.set_linear_coefficients(
        {rxn.forward_variable: 1, rxn.reverse_variable: -1})
    _model.slim_optimize()
    sutil.check_solver_status(_model.solver.status)
    if _loopless:
        value = loopless_fva_iter(_model, rxn)
    else:
        value = _model.solver.objective.value
    _model.solver.objective.set_linear_coefficients(
        {rxn.forward_variable: 0, rxn.reverse_variable: 0})
    return reaction_id, value


def flux_variability_analysis(model, reaction_list=None, loopless=False,
                              fraction_of_optimum=1.0, pfba_factor=None,
                              processes=None):
    """
    Determine the minimum and maximum possible flux value for each reaction.

    Parameters
    ----------
    model : cobra.Model
        The model for which to run the analysis. It will *not* be modified.
    reaction_list : list of cobra.Reaction or str, optional
        The reactions for which to obtain min/max fluxes. If None will use
        all reactions in the model (default).
    loopless : boolean, optional
        Whether to return only loopless solutions. This is significantly
        slower. Please also refer to the notes.
    fraction_of_optimum : float, optional
        Must be <= 1.0. Requires that the objective value is at least the
        fraction times maximum objective value. A value of 0.85 for instance
        means that the objective has to be at least at 85% percent of its
        maximum.
    pfba_factor : float, optional
        Add an additional constraint to the model that requires the total sum
        of absolute fluxes must not be larger than this value times the
        smallest possible sum of absolute fluxes, i.e., by setting the value
        to 1.1 the total sum of absolute fluxes must not be more than
        10% larger than the pFBA solution. Since the pFBA solution is the
        one that optimally minimizes the total flux sum, the ``pfba_factor``
        should, if set, be larger than one. Setting this value may lead to
        more realistic predictions of the effective flux bounds.
    processes : int, optional
        The number of parallel processes to run. If not explicitly passed,
        will be set from the global configuration singleton.

    Returns
    -------
    pandas.DataFrame
        A data frame with reaction identifiers as the index and two columns:
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
    "naive" version using ``add_loopless(model)``. Also note that if you have
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
        reaction_ids = [r.id for r in model.reactions]
    else:
        reaction_ids = [r.id
                        for r in model.reactions.get_by_any(reaction_list)]

    if processes is None:
        processes = CONFIGURATION.processes
    num_reactions = len(reaction_ids)
    processes = min(processes, num_reactions)

    fva_result = DataFrame({
        "minimum": zeros(num_reactions, dtype=float),
        "maximum": zeros(num_reactions, dtype=float)
    }, index=reaction_ids)
    prob = model.problem
    with model:
        # Safety check before setting up FVA.
        model.slim_optimize(error_value=None,
                            message="There is no optimal solution for the "
                                    "chosen objective!")
        # Add the previous objective as a variable to the model then set it to
        # zero. This also uses the fraction to create the lower/upper bound for
        # the old objective.
        # TODO: Use utility function here (fix_objective_as_constraint)?
        if model.solver.objective.direction == "max":
            fva_old_objective = prob.Variable(
                "fva_old_objective",
                lb=fraction_of_optimum * model.solver.objective.value)
        else:
            fva_old_objective = prob.Variable(
                "fva_old_objective",
                ub=fraction_of_optimum * model.solver.objective.value)
        fva_old_obj_constraint = prob.Constraint(
            model.solver.objective.expression - fva_old_objective, lb=0, ub=0,
            name="fva_old_objective_constraint")
        model.add_cons_vars([fva_old_objective, fva_old_obj_constraint])

        if pfba_factor is not None:
            if pfba_factor < 1.:
                warn("The 'pfba_factor' should be larger or equal to 1.",
                     UserWarning)
            with model:
                add_pfba(model, fraction_of_optimum=0)
                ub = model.slim_optimize(error_value=None)
                flux_sum = prob.Variable("flux_sum", ub=pfba_factor * ub)
                flux_sum_constraint = prob.Constraint(
                    model.solver.objective.expression - flux_sum, lb=0, ub=0,
                    name="flux_sum_constraint")
            model.add_cons_vars([flux_sum, flux_sum_constraint])

        model.objective = Zero  # This will trigger the reset as well
        for what in ("minimum", "maximum"):
            if processes > 1:
                # We create and destroy a new pool here in order to set the
                # objective direction for all reactions. This creates a
                # slight overhead but seems the most clean.
                chunk_size = len(reaction_ids) // processes
                pool = multiprocessing.Pool(
                    processes,
                    initializer=_init_worker,
                    initargs=(model, loopless, what[:3])
                )
                for rxn_id, value in pool.imap_unordered(_fva_step,
                                                         reaction_ids,
                                                         chunksize=chunk_size):
                    fva_result.at[rxn_id, what] = value
                pool.close()
                pool.join()
            else:
                _init_worker(model, loopless, what[:3])
                for rxn_id, value in map(_fva_step, reaction_ids):
                    fva_result.at[rxn_id, what] = value

    return fva_result[["minimum", "maximum"]]


def find_blocked_reactions(model,
                           reaction_list=None,
                           zero_cutoff=None,
                           open_exchanges=False,
                           processes=None):
    """
    Find reactions that cannot carry any flux.

    The question whether or not a reaction is blocked is highly dependent
    on the current exchange reaction settings for a COBRA model. Hence an
    argument is provided to open all exchange reactions.

    Notes
    -----
    Sink and demand reactions are left untouched. Please modify them manually.

    Parameters
    ----------
    model : cobra.Model
        The model to analyze.
    reaction_list : list, optional
        List of reactions to consider, the default includes all model
        reactions.
    zero_cutoff : float, optional
        Flux value which is considered to effectively be zero
        (default model.tolerance).
    open_exchanges : bool, optional
        Whether or not to open all exchange reactions to very high flux ranges.
    processes : int, optional
        The number of parallel processes to run. Can speed up the computations
        if the number of reactions is large. If not explicitly
        passed, it will be set from the global configuration singleton.

    Returns
    -------
    list
        List with the identifiers of blocked reactions.

    """
    zero_cutoff = normalize_cutoff(model, zero_cutoff)

    with model:
        if open_exchanges:
            for reaction in model.exchanges:
                reaction.bounds = (min(reaction.lower_bound, -1000),
                                   max(reaction.upper_bound, 1000))
        if reaction_list is None:
            reaction_list = model.reactions
        # Limit the search space to reactions which have zero flux. If the
        # reactions already carry flux in this solution,
        # then they cannot be blocked.
        model.slim_optimize()
        solution = get_solution(model, reactions=reaction_list)
        reaction_list = solution.fluxes[
            solution.fluxes.abs() < zero_cutoff].index.tolist()
        # Run FVA to find reactions where both the minimal and maximal flux
        # are zero (below the cut off).
        flux_span = flux_variability_analysis(
            model, fraction_of_optimum=0., reaction_list=reaction_list,
            processes=processes
        )
        return flux_span[
            flux_span.abs().max(axis=1) < zero_cutoff].index.tolist()


def find_essential_genes(model, threshold=None, processes=None):
    """
    Return a set of essential genes.

    A gene is considered essential if restricting the flux of all reactions
    that depend on it to zero causes the objective, e.g., the growth rate,
    to also be zero, below the threshold, or infeasible.

    Parameters
    ----------
    model : cobra.Model
        The model to find the essential genes for.
    threshold : float, optional
        Minimal objective flux to be considered viable. By default this is
        1% of the maximal objective.
    processes : int, optional
        The number of parallel processes to run. If not passed,
        will be set to the number of CPUs found.
    processes : int, optional
        The number of parallel processes to run. Can speed up the computations
        if the number of knockouts to perform is large. If not explicitly
        passed, it will be set from the global configuration singleton.

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
    return {model.genes.get_by_id(g) for ids in essential for g in ids}


def find_essential_reactions(model, threshold=None, processes=None):
    """Return a set of essential reactions.

    A reaction is considered essential if restricting its flux to zero
    causes the objective, e.g., the growth rate, to also be zero, below the
    threshold, or infeasible.


    Parameters
    ----------
    model : cobra.Model
        The model to find the essential reactions for.
    threshold : float, optional
        Minimal objective flux to be considered viable. By default this is
        1% of the maximal objective.
    processes : int, optional
        The number of parallel processes to run. Can speed up the computations
        if the number of knockouts to perform is large. If not explicitly
        passed, it will be set from the global configuration singleton.

    Returns
    -------
    set
        Set of essential reactions
    """
    if threshold is None:
        threshold = model.slim_optimize(error_value=None) * 1E-02
    deletions = single_reaction_deletion(
        model, method='fba', processes=processes)
    essential = deletions.loc[deletions['growth'].isna() |
                              (deletions['growth'] < threshold), :].index
    return {model.reactions.get_by_id(r) for ids in essential for r in ids}
