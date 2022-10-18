"""Provide variability based methods such as flux variability or gene essentiality."""


import logging
from typing import TYPE_CHECKING, List, Optional, Set, Tuple, Union
from warnings import warn

import numpy as np
import pandas as pd
from optlang.symbolics import Zero

from ..core import Configuration, get_solution
from ..util import ProcessPool
from ..util import solver as sutil
from .deletion import single_gene_deletion, single_reaction_deletion
from .helpers import normalize_cutoff
from .loopless import loopless_fva_iter
from .parsimonious import add_pfba


if TYPE_CHECKING:
    from cobra import Gene, Model, Reaction


logger = logging.getLogger(__name__)
configuration = Configuration()


def _init_worker(model: "Model", loopless: bool, sense: str) -> None:
    """Initialize a global model object for multiprocessing.

    Parameters
    ----------
    model: cobra.Model
        The model to operate on.
    loopless: bool
        Whether to use loopless version.
    sense: {"max", "min"}
        Whether to maximise or minimise objective.

    """
    global _model
    global _loopless
    _model = model
    _model.solver.objective.direction = sense
    _loopless = loopless


def _fva_step(reaction_id: str) -> Tuple[str, float]:
    """Take a step for calculating FVA.

    Parameters
    ----------
    reaction_id: str
        The ID of the reaction.

    Returns
    -------
    tuple of (str, float)
        The reaction ID with the flux value.

    """
    global _model
    global _loopless
    rxn = _model.reactions.get_by_id(reaction_id)
    # The previous objective assignment already triggers a reset
    # so directly update coefs here to not trigger redundant resets
    # in the history manager which can take longer than the actual
    # FVA for small models
    _model.solver.objective.set_linear_coefficients(
        {rxn.forward_variable: 1, rxn.reverse_variable: -1}
    )
    _model.slim_optimize()
    sutil.check_solver_status(_model.solver.status)
    if _loopless:
        value = loopless_fva_iter(_model, rxn)
    else:
        value = _model.solver.objective.value
    # handle infeasible case
    if value is None:
        value = float("nan")
        logger.warning(
            f"Could not get flux for reaction {rxn.id}, setting it to NaN. "
            "This is usually due to numerical instability."
        )
    _model.solver.objective.set_linear_coefficients(
        {rxn.forward_variable: 0, rxn.reverse_variable: 0}
    )
    return reaction_id, value


def flux_variability_analysis(
    model: "Model",
    reaction_list: Optional[List[Union["Reaction", str]]] = None,
    loopless: bool = False,
    fraction_of_optimum: float = 1.0,
    pfba_factor: Optional[float] = None,
    processes: Optional[int] = None,
) -> pd.DataFrame:
    """Determine the minimum and maximum flux value for each reaction.

    Parameters
    ----------
    model : cobra.Model
        The model for which to run the analysis. It will *not* be modified.
    reaction_list : list of cobra.Reaction or str, optional
        The reactions for which to obtain min/max fluxes. If None will use
        all reactions in the model (default None).
    loopless : bool, optional
        Whether to return only loopless solutions. This is significantly
        slower. Please also refer to the notes (default False).
    fraction_of_optimum : float, optional
        Must be <= 1.0. Requires that the objective value is at least the
        fraction times maximum objective value. A value of 0.85 for instance
        means that the objective has to be at least at 85% percent of its
        maximum (default 1.0).
    pfba_factor : float, optional
        Add an additional constraint to the model that requires the total sum
        of absolute fluxes must not be larger than this value times the
        smallest possible sum of absolute fluxes, i.e., by setting the value
        to 1.1 the total sum of absolute fluxes must not be more than
        10% larger than the pFBA solution. Since the pFBA solution is the
        one that optimally minimizes the total flux sum, the `pfba_factor`
        should, if set, be larger than one. Setting this value may lead to
        more realistic predictions of the effective flux bounds
        (default None).
    processes : int, optional
        The number of parallel processes to run. If not explicitly passed,
        will be set from the global configuration singleton (default None).

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
    sub-optimal.

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
        reaction_ids = [r.id for r in model.reactions]
    else:
        reaction_ids = [r.id for r in model.reactions.get_by_any(reaction_list)]

    if processes is None:
        processes = configuration.processes

    num_reactions = len(reaction_ids)
    processes = min(processes, num_reactions)

    fva_result = pd.DataFrame(
        {
            "minimum": np.zeros(num_reactions, dtype=float),
            "maximum": np.zeros(num_reactions, dtype=float),
        },
        index=reaction_ids,
    )
    prob = model.problem
    with model:
        # Safety check before setting up FVA.
        model.slim_optimize(
            error_value=None,
            message="There is no optimal solution for the chosen objective!",
        )
        # Add the previous objective as a variable to the model then set it to
        # zero. This also uses the fraction to create the lower/upper bound for
        # the old objective.
        # TODO: Use utility function here (fix_objective_as_constraint)?
        if model.solver.objective.direction == "max":
            fva_old_objective = prob.Variable(
                "fva_old_objective",
                lb=fraction_of_optimum * model.solver.objective.value,
            )
        else:
            fva_old_objective = prob.Variable(
                "fva_old_objective",
                ub=fraction_of_optimum * model.solver.objective.value,
            )
        fva_old_obj_constraint = prob.Constraint(
            model.solver.objective.expression - fva_old_objective,
            lb=0,
            ub=0,
            name="fva_old_objective_constraint",
        )
        model.add_cons_vars([fva_old_objective, fva_old_obj_constraint])

        if pfba_factor is not None:
            if pfba_factor < 1.0:
                warn(
                    "The 'pfba_factor' should be larger or equal to 1.",
                    UserWarning,
                )
            with model:
                add_pfba(model, fraction_of_optimum=0)
                ub = model.slim_optimize(error_value=None)
                flux_sum = prob.Variable("flux_sum", ub=pfba_factor * ub)
                flux_sum_constraint = prob.Constraint(
                    model.solver.objective.expression - flux_sum,
                    lb=0,
                    ub=0,
                    name="flux_sum_constraint",
                )
            model.add_cons_vars([flux_sum, flux_sum_constraint])

        model.objective = Zero  # This will trigger the reset as well
        for what in ("minimum", "maximum"):
            if processes > 1:
                # We create and destroy a new pool here in order to set the
                # objective direction for all reactions. This creates a
                # slight overhead but seems the most clean.
                chunk_size = len(reaction_ids) // processes
                with ProcessPool(
                    processes,
                    initializer=_init_worker,
                    initargs=(model, loopless, what[:3]),
                ) as pool:
                    for rxn_id, value in pool.imap_unordered(
                        _fva_step, reaction_ids, chunksize=chunk_size
                    ):
                        fva_result.at[rxn_id, what] = value
            else:
                _init_worker(model, loopless, what[:3])
                for rxn_id, value in map(_fva_step, reaction_ids):
                    fva_result.at[rxn_id, what] = value

    return fva_result[["minimum", "maximum"]]


def find_blocked_reactions(
    model: "Model",
    reaction_list: Optional[List[Union["Reaction", str]]] = None,
    zero_cutoff: Optional[float] = None,
    open_exchanges: bool = False,
    processes: Optional[int] = None,
) -> List["Reaction"]:
    """Find reactions that cannot carry any flux.

    The question whether or not a reaction is blocked is highly dependent
    on the current exchange reaction settings for a COBRA model. Hence an
    argument is provided to open all exchange reactions.

    Parameters
    ----------
    model : cobra.Model
        The model to analyze.
    reaction_list : list of cobra.Reaction or str, optional
        List of reactions to consider, the default includes all model
        reactions (default None).
    zero_cutoff : float, optional
        Flux value which is considered to effectively be zero. The default
        is set to use `model.tolerance` (default None).
    open_exchanges : bool, optional
        Whether or not to open all exchange reactions to very high flux
        ranges (default False).
    processes : int, optional
        The number of parallel processes to run. Can speed up the
        computations if the number of reactions is large. If not explicitly
        passed, it will be set from the global configuration singleton
        (default None).

    Returns
    -------
    list of cobra.Reaction
        List with the identifiers of blocked reactions.

    Notes
    -----
    Sink and demand reactions are left untouched. Please modify them manually.

    """
    zero_cutoff = normalize_cutoff(model, zero_cutoff)

    with model:
        if open_exchanges:
            for reaction in model.exchanges:
                reaction.bounds = (
                    min(reaction.lower_bound, -1000),
                    max(reaction.upper_bound, 1000),
                )
        if reaction_list is None:
            reaction_list = model.reactions
        # Limit the search space to reactions which have zero flux. If the
        # reactions already carry flux in this solution,
        # then they cannot be blocked.
        model.slim_optimize()
        solution = get_solution(model, reactions=reaction_list)
        reaction_list = solution.fluxes[
            solution.fluxes.abs() < zero_cutoff
        ].index.tolist()
        # Run FVA to find reactions where both the minimal and maximal flux
        # are zero (below the cut off).
        flux_span = flux_variability_analysis(
            model,
            fraction_of_optimum=0.0,
            reaction_list=reaction_list,
            processes=processes,
        )
        return flux_span[flux_span.abs().max(axis=1) < zero_cutoff].index.tolist()


def find_essential_genes(
    model: "Model",
    threshold: Optional[float] = None,
    processes: Optional[int] = None,
) -> Set["Gene"]:
    """Return a set of essential genes.

    A gene is considered essential if restricting the flux of all reactions
    that depend on it to zero causes the objective, e.g., the growth rate,
    to also be zero, below the threshold, or infeasible.

    Parameters
    ----------
    model : cobra.Model
        The model to find the essential genes for.
    threshold : float, optional
        Minimal objective flux to be considered viable. By default this is
        1% of the maximal objective (default None).
    processes : int, optional
        The number of parallel processes to run. Can speed up the computations
        if the number of knockouts to perform is large. If not explicitly
        passed, it will be set from the global configuration singleton
        (default None).

    Returns
    -------
    set of cobra.Gene
        Set of essential genes.

    """
    if threshold is None:
        threshold = model.slim_optimize(error_value=None) * 1e-02
    deletions = single_gene_deletion(model, method="fba", processes=processes)
    essential = deletions.loc[
        deletions["growth"].isna() | (deletions["growth"] < threshold), :
    ].ids
    return {model.genes.get_by_id(g) for ids in essential for g in ids}


def find_essential_reactions(
    model: "Model",
    threshold: Optional[float] = None,
    processes: Optional[int] = None,
) -> Set["Reaction"]:
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
        1% of the maximal objective (default None).
    processes : int, optional
        The number of parallel processes to run. Can speed up the computations
        if the number of knockouts to perform is large. If not explicitly
        passed, it will be set from the global configuration singleton
        (default None).

    Returns
    -------
    set of cobra.Reaction
        Set of essential reactions.

    """
    if threshold is None:
        threshold = model.slim_optimize(error_value=None) * 1e-02
    deletions = single_reaction_deletion(model, method="fba", processes=processes)
    essential = deletions.loc[
        deletions["growth"].isna() | (deletions["growth"] < threshold), :
    ].ids
    return {model.reactions.get_by_id(r) for ids in essential for r in ids}


def faster_fva(
    model: "Model",
    reaction_list: Optional[List[Union["Reaction", str]]] = None,
    loopless: bool = False,
    fraction_of_optimum: float = 1.0,
    pfba_factor: Optional[float] = None,
    tol: Optional[float] = 1e-13,
) -> pd.DataFrame:
    """Determine the minimum and maximum flux value for each reaction.

    This function uses a faster algorithm than flux_variability_analysis.

    Parameters
    ----------
    model : cobra.Model
        The model for which to run the analysis. It will *not* be modified.
    reaction_list : list of cobra.Reaction or str, optional
        The reactions for which to obtain min/max fluxes. If None will use
        all reactions in the model (default None).
    loopless : bool, optional
        Whether to return only loopless solutions. This is significantly
        slower. Please also refer to the notes (default False).
    fraction_of_optimum : float, optional
        Must be <= 1.0. Requires that the objective value is at least the
        fraction times maximum objective value. A value of 0.85 for instance
        means that the objective has to be at least at 85% percent of its
        maximum (default 1.0).
    pfba_factor : float, optional
        Add an additional constraint to the model that requires the total sum
        of absolute fluxes must not be larger than this value times the
        smallest possible sum of absolute fluxes, i.e., by setting the value
        to 1.1 the total sum of absolute fluxes must not be more than
        10% larger than the pFBA solution. Since the pFBA solution is the
        one that optimally minimizes the total flux sum, the `pfba_factor`
        should, if set, be larger than one. Setting this value may lead to
        more realistic predictions of the effective flux bounds
        (default None).
    tol : float, optional.
        Tolerance to consider that a flux value is equal to the lower or
        upper bounds (default 1e-13).

    Returns
    -------
    pandas.DataFrame
        A data frame with reaction identifiers as the index and two columns:
        - maximum: indicating the highest possible flux
        - minimum: indicating the lowest possible flux

    Notes
    -----
    This implements the a faster algorithm described in [1]_. Please note that
    the flux distribution containing all minimal/maximal fluxes does not have
    to be a feasible solution for the model. Fluxes are minimized/maximized
    individually and a single minimal flux might require all others to be
    sub-optimal.

    Using the loopless option will lead to a significant increase in
    computation time (about a factor of 100 for large models). However, the
    algorithm used here (see [2]_) is still more than 1000x faster than the
    "naive" version using `add_loopless(model)`. Also note that if you have
    included constraints that force a loop (for instance by setting all fluxes
    in a loop to be non-zero) this loop will be included in the solution.

    References
    ----------
    .. [1] An Improved Algorithm for Flux Variability Analysis
       Kenefake D, Armingol E, Lewis NE, Pistikopoulos EN.
       Research Square (Pre-Print)
       doi: https://doi.org/10.21203/rs.3.rs-2109496/v1

    .. [2] CycleFreeFlux: efficient removal of thermodynamically infeasible
       loops from flux distributions.
       Desouki AA, Jarre F, Gelius-Dietrich G, Lercher MJ.
       Bioinformatics. 2015 Jul 1;31(13):2159-65.
       doi: 10.1093/bioinformatics/btv096.

    """
    if reaction_list is None:
        reaction_ids = [r.id for r in model.reactions]
    else:
        reaction_ids = [r.id for r in model.reactions.get_by_any(reaction_list)]

    num_reactions = len(reaction_ids)

    fva_result = pd.DataFrame(
        {
            "minimum": np.empty(num_reactions, dtype=float).fill(np.nan),
            "maximum": np.empty(num_reactions, dtype=float).fill(np.nan),
        },
        index=reaction_ids,
        dtype="float64",
    )

    lower_bound_problems = set(range(len(reaction_ids)))
    upper_bound_problems = set(range(len(reaction_ids)))
    prob = model.problem

    def remove_bound_problems(v_vals, lower_bounds, upper_bounds):
        """Remove problems to solve from the lower bound.

        This is helper function.
        """
        lbp_to_remove = np.where(
            (v_vals <= lower_bounds + tol) & (v_vals >= lower_bounds - tol)
        )[0]
        lower_bound_problems.difference_update(set(lbp_to_remove))
        if len(lbp_to_remove) > 0:
            fva_result.loc[
                np.array(reaction_ids)[lbp_to_remove], "minimum"
            ] = lower_bounds[lbp_to_remove]

        ubp_to_remove = np.where(
            (v_vals <= upper_bounds + tol) & (v_vals >= upper_bounds - tol)
        )[0]
        upper_bound_problems.difference_update(set(ubp_to_remove))
        if len(ubp_to_remove) > 0:
            fva_result.loc[
                np.array(reaction_ids)[ubp_to_remove], "maximum"
            ] = upper_bounds[ubp_to_remove]

    with model:
        # Safety check before setting up FVA.
        model.slim_optimize(
            error_value=None,
            message="There is no optimal solution for the chosen objective!",
        )

        # Find reactions that hit bounds when optimizing original objective
        v_vals = np.array([model.reactions.get_by_id(idx).flux for idx in reaction_ids])
        lower_bounds = np.array(
            [model.reactions.get_by_id(idx).lower_bound for idx in reaction_ids]
        )
        upper_bounds = np.array(
            [model.reactions.get_by_id(idx).upper_bound for idx in reaction_ids]
        )
        remove_bound_problems(v_vals, lower_bounds, upper_bounds)

        number_lps = 1

        # Add the previous objective as a variable to the model then set it to
        # zero. This also uses the fraction to create the lower/upper bound for
        # the old objective.
        # TODO: Use utility function here (fix_objective_as_constraint)?
        if model.solver.objective.direction == "max":
            fva_old_objective = prob.Variable(
                "fva_old_objective",
                lb=fraction_of_optimum * model.solver.objective.value,
            )
        else:
            fva_old_objective = prob.Variable(
                "fva_old_objective",
                ub=fraction_of_optimum * model.solver.objective.value,
            )
        fva_old_obj_constraint = prob.Constraint(
            model.solver.objective.expression - fva_old_objective,
            lb=0,
            ub=0,
            name="fva_old_objective_constraint",
        )
        model.add_cons_vars([fva_old_objective, fva_old_obj_constraint])

        if pfba_factor is not None:
            if pfba_factor < 1.0:
                warn(
                    "The 'pfba_factor' should be larger or equal to 1.",
                    UserWarning,
                )
            with model:
                add_pfba(model, fraction_of_optimum=0)
                ub = model.slim_optimize(error_value=None)
                flux_sum = prob.Variable("flux_sum", ub=pfba_factor * ub)
                flux_sum_constraint = prob.Constraint(
                    model.solver.objective.expression - flux_sum,
                    lb=0,
                    ub=0,
                    name="flux_sum_constraint",
                )
            model.add_cons_vars([flux_sum, flux_sum_constraint])

        model.objective = Zero  # This will trigger the reset as well

        # solve the lower bound problems
        what = "minimum"
        while len(lower_bound_problems) != 0:
            # Grab a problem of the top
            lb_problem = lower_bound_problems.pop()
            rxn_id = reaction_ids[lb_problem]

            # Run FVA step
            _init_worker(model, loopless, what[:3])
            value = _fva_step(rxn_id)
            fva_result.at[rxn_id, what] = value[1]
            number_lps += 1

            # Find reactions that hit bounds when optimizing objective
            v_vals = np.array(
                [model.reactions.get_by_id(idx).flux for idx in reaction_ids]
            )
            lower_bounds = np.array(
                [model.reactions.get_by_id(idx).lower_bound for idx in reaction_ids]
            )
            upper_bounds = np.array(
                [model.reactions.get_by_id(idx).upper_bound for idx in reaction_ids]
            )
            remove_bound_problems(v_vals, lower_bounds, upper_bounds)

        # solve the upper bound problems
        what = "maximum"
        while len(upper_bound_problems) != 0:
            # grab a problem of the top
            ub_problem = upper_bound_problems.pop()
            rxn_id = reaction_ids[ub_problem]

            # Run FVA step
            _init_worker(model, loopless, what[:3])
            value = _fva_step(rxn_id)
            fva_result.at[rxn_id, what] = value[1]
            number_lps += 1

            # Find reactions that hit bounds when optimizing objective
            v_vals = np.array(
                [model.reactions.get_by_id(idx).flux for idx in reaction_ids]
            )
            lower_bounds = np.array(
                [model.reactions.get_by_id(idx).lower_bound for idx in reaction_ids]
            )
            upper_bounds = np.array(
                [model.reactions.get_by_id(idx).upper_bound for idx in reaction_ids]
            )
            remove_bound_problems(v_vals, lower_bounds, upper_bounds)
    # total_lps = 2 * len(reaction_ids) + 1
    # print(
    #     f"{number_lps} LPs were solved during this analysis. \
    # The regular algorithm would have solved {total_lps} LPs."
    # )
    return fva_result[["minimum", "maximum"]]
