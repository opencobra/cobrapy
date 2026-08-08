"""Provide variability based methods such as flux variability or gene essentiality."""

import logging
from typing import TYPE_CHECKING, Dict, List, Optional, Set, Tuple, Union
from warnings import warn

import numpy as np
import pandas as pd
from optlang.interface import OPTIMAL, TIME_LIMIT
from optlang.symbolics import Zero

from ..core import Configuration, get_solution
from ..util import ProcessPool
from ..util import solver as sutil
from .deletion import single_gene_deletion, single_reaction_deletion
from .find_cyclic_reactions import find_cyclic_reactions
from .helpers import normalize_cutoff
from .loopless import add_loopless, loopless_fva_iter
from .parsimonious import add_pfba


if TYPE_CHECKING:
    from cobra import Gene, Model, Reaction


logger = logging.getLogger(__name__)
configuration = Configuration()


def _init_worker(
    model: "Model",
    loopless: bool,
    sense: str,
    return_fluxes: bool = False,
    time_limit: Optional[float] = None,
    accept_incumbent: bool = False,
) -> None:
    """Initialize a global model object for multiprocessing.

    Parameters
    ----------
    model: cobra.Model
        The model to operate on.
    loopless: bool
        Whether to use loopless version.
    sense: {"max", "min"}
        Whether to maximise or minimise objective.
    return_fluxes: bool
        Whether to keep the flux distribution of each optimum.
    time_limit: float, optional
        Per-solve time limit in seconds, applied to this worker's solver
        (default None, no limit).
    accept_incumbent: bool
        Whether a solve that ends at the time limit records the solver's
        best feasible solution instead of NaN.

    """
    global _model
    global _loopless
    global _return_fluxes
    global _accept_incumbent
    _model = model
    _model.solver.objective.direction = sense
    _loopless = loopless
    _return_fluxes = return_fluxes
    _accept_incumbent = accept_incumbent
    if time_limit is not None:
        model.solver.configuration.timeout = time_limit


def _fva_step(
    reaction_id: str,
) -> Tuple[str, float, Optional[Dict[str, float]], bool]:
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
    global _return_fluxes
    global _accept_incumbent
    rxn = _model.reactions.get_by_id(reaction_id)
    # The previous objective assignment already triggers a reset
    # so directly update coefs here to not trigger redundant resets
    # in the history manager which can take longer than the actual
    # FVA for small models
    _model.solver.objective.set_linear_coefficients(
        {rxn.forward_variable: 1, rxn.reverse_variable: -1}
    )
    _model.slim_optimize()
    status = _model.solver.status
    optimal = status == OPTIMAL
    fluxes = None
    value = None
    if _loopless:
        if not optimal:
            sutil.check_solver_status(status)
        elif _return_fluxes:
            solution = loopless_fva_iter(_model, rxn, solution=True)
            value = None if solution is None else solution.fluxes[reaction_id]
            fluxes = None if solution is None else solution.fluxes
        else:
            value = loopless_fva_iter(_model, rxn)
    elif optimal:
        value = _model.solver.objective.value
        if _return_fluxes:
            fluxes = get_solution(_model).fluxes
    elif _accept_incumbent and status == TIME_LIMIT:
        # The solver holds its best feasible solution (the incumbent) when a
        # time limit interrupts branch and bound. Record it, flagged as
        # unproven; if no incumbent exists yet, extraction raises and the
        # value stays NaN.
        try:
            value = _model.solver.objective.value
            if _return_fluxes:
                fluxes = get_solution(_model).fluxes
        except Exception:
            value = None
            fluxes = None
    else:
        sutil.check_solver_status(status)
    # handle infeasible case
    if value is None or value != value:
        value = float("nan")
        fluxes = None
        logger.warning(
            f"Could not get flux for reaction {rxn.id}, setting it to NaN. "
            "This is usually due to numerical instability or a time limit "
            "hit before any feasible solution was found."
        )
    _model.solver.objective.set_linear_coefficients(
        {rxn.forward_variable: 0, rxn.reverse_variable: 0}
    )
    proven = optimal and value == value
    return reaction_id, value, fluxes, proven


def flux_variability_analysis(
    model: "Model",
    reaction_list: Optional[List[Union["Reaction", str]]] = None,
    loopless: Union[Optional[str], bool] = None,
    fraction_of_optimum: float = 1.0,
    pfba_factor: Optional[float] = None,
    processes: Optional[int] = None,
    return_fluxes: bool = False,
    time_limit: Optional[float] = None,
    accept_incumbent: bool = False,
) -> Union[pd.DataFrame, Tuple[pd.DataFrame, Dict[str, pd.DataFrame]]]:
    """Determine the minimum and maximum flux value for each reaction.

    Parameters
    ----------
    model : cobra.Model
        The model for which to run the analysis. It will *not* be modified.
    reaction_list : list of cobra.Reaction or str, optional
        The reactions for which to obtain min/max fluxes. If None will use
        all reactions in the model (default None).
    loopless : str, "fastSNP" or "cycleFreeFlux", optional
        If this value is set, only loopless solutions will be returned.
        Boolean values are deprecated. Provided value means the algorithm
        to constrain the model to loopless solutions.
        Please also refer to the notes (default None).
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
    return_fluxes : bool, optional
        Whether to also return the flux distribution found at each optimum
        (default False). FVA solves one optimization per reaction per direction
        and discards every solution but the objective value; setting this keeps
        them, which saves recomputing an FVA-sized batch of solves when the
        distributions themselves are wanted -- for example as a starting pool
        for sampling. With `loopless="fastSNP"` the loop constraints are applied
        to the whole model when this is set, so that every returned distribution
        is loopless; this is slower than the default, which constrains only the
        reactions that can carry a loop.
    time_limit : float, optional
        Per-solve time limit in seconds, applied to every optimization in the
        analysis (default None, no limit). Mainly useful with the MILP
        loopless variants ("fastSNP", "potentials") at genome scale, where
        individual solves may not prove optimality in reasonable time. What a
        solve that hits the limit records depends on `accept_incumbent`.
    accept_incumbent : bool, optional
        Only meaningful for solves that end non-optimal at the `time_limit`.
        If False (default), such solves record NaN. If True, the solver's
        best feasible solution found so far (the incumbent) is recorded
        instead, and the result gains boolean columns ``minimum_proven`` /
        ``maximum_proven`` marking, per direction, whether the value carries
        an optimality proof. Unproven values are one-sided estimates: an
        unproven maximum is a lower bound on the true maximum, and vice
        versa. With `return_fluxes`, incumbent flux distributions are
        returned like optimal ones -- they are feasible (loopless under the
        MILP variants, up to the solver's integrality tolerance) but not
        necessarily extreme. Ignored on the "cycleFreeFlux" path.

    Returns
    -------
    pandas.DataFrame or tuple of (pandas.DataFrame, dict of {str: pandas.DataFrame})
        A data frame with reaction identifiers as the index and two columns:
        - maximum: indicating the highest possible flux
        - minimum: indicating the lowest possible flux
        If `time_limit` is set, two additional boolean columns
        ``minimum_proven`` and ``maximum_proven`` report per direction
        whether the recorded value was solved to proven optimality.
        If `return_fluxes` is True, a tuple whose second element maps
        "minimum" and "maximum" to data frames of the flux distributions, each
        indexed by the optimized reaction with reaction identifiers as columns.

    Notes
    -----
    This implements the fast version as described in [1]_. Please note that
    the flux distribution containing all minimal/maximal fluxes does not have
    to be a feasible solution for the model. Fluxes are minimized/maximized
    individually and a single minimal flux might require all others to be
    sub-optimal.

    Using the loopless option will lead to a significant increase in
    computation time (about a factor of 100 for large models).

    If `loopless` is set to "fastSNP", the optimal loopless flux bounds will be
    found by adding the loopless constraints to the model using efficient
    Fast-SNP algorithm (see [2]_).

    If `loopless` is set to "cycleFreeFlux", the loops removal algorithm will be
    used (see [3]_). Note: this algorithm does not guarantee to find optimal bounds.

    References
    ----------
    .. [1] Computationally efficient flux variability analysis.
       Gudmundsson S, Thiele I.
       BMC Bioinformatics. 2010 Sep 29;11:489.
       doi: 10.1186/1471-2105-11-489, PMID: 20920235

    .. [2] Fast-SNP: a fast matrix pre-processing algorithm for efficient
       loopless flux optimization of metabolic models. Saa PA, Nielsen LK.
       Bioinformatics. 2016 Dec;32(24):3807–3814. doi: 10.1093/bioinformatics/btw555.

    .. [3] CycleFreeFlux: efficient removal of thermodynamically infeasible
       loops from flux distributions.
       Desouki AA, Jarre F, Gelius-Dietrich G, Lercher MJ.
       Bioinformatics. 2015 Jul 1;31(13):2159-65.
       doi: 10.1093/bioinformatics/btv096.
    """
    if loopless is not None and isinstance(loopless, bool):
        warn(
            "Passing a boolean value to the `loopless` argument is deprecated. "
            "Please pass either None, 'fastSNP' or 'cycleFreeFlux'.",
            DeprecationWarning,
            stacklevel=2,
        )
        loopless = "cycleFreeFlux" if loopless else None

    if loopless not in (None, "fastSNP", "cycleFreeFlux"):
        raise ValueError(
            "The `loopless` argument must be either None, 'fastSNP' or 'cycleFreeFlux'."
        )

    if reaction_list is None:
        reaction_ids = [r.id for r in model.reactions]
    else:
        reaction_ids = [r.id for r in model.reactions.get_by_any(reaction_list)]

    if processes is None:
        processes = configuration.processes

    num_reactions = len(reaction_ids)
    processes = min(processes, num_reactions)
    # In serial mode _init_worker runs on the caller's model; remember the
    # solver timeout so it can be restored afterwards.
    orig_timeout = model.solver.configuration.timeout

    fva_result = pd.DataFrame(
        {
            "minimum": np.zeros(num_reactions, dtype=float),
            "maximum": np.zeros(num_reactions, dtype=float),
            "minimum_proven": np.ones(num_reactions, dtype=bool),
            "maximum_proven": np.ones(num_reactions, dtype=bool),
        },
        index=reaction_ids,
    )
    result_columns = ["minimum", "maximum"]
    if time_limit is not None:
        result_columns += ["minimum_proven", "maximum_proven"]

    reaction_ids_by_type = [
        {
            "minimum": [],
            "maximum": [],
        },
        {
            "minimum": [],
            "maximum": [],
        },
    ]
    if loopless is not None:
        cyclic_reactions, cyclic_directions = find_cyclic_reactions(model)
        cyclic_reaction_index = {r_id: i for i, r_id in enumerate(cyclic_reactions)}
        for r_id in reaction_ids:
            i = cyclic_reaction_index.get(r_id)
            for loc, dir in enumerate(("minimum", "maximum")):
                if i is not None and cyclic_directions[i][loc]:
                    reaction_ids_by_type[1][dir].append(r_id)
                else:
                    reaction_ids_by_type[0][dir].append(r_id)
    else:
        reaction_ids_by_type[0]["minimum"] = reaction_ids
        reaction_ids_by_type[0]["maximum"] = reaction_ids

    flux_rows: Dict[str, Dict[str, pd.Series]] = {"minimum": {}, "maximum": {}}

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

        # Reactions that cannot carry a loop are normally optimized without the
        # loop constraints: their optimum is the same either way, so only the
        # objective value is needed and the LP is cheaper. But that shortcut also
        # means their solution VECTORS may contain loops elsewhere in the network.
        # When the caller asks for those vectors, constrain the whole model up
        # front so every returned distribution is loopless. The bounds are
        # unaffected -- the loop law over cyclic reactions is the entire loop law,
        # so the feasible set is identical either way.
        constrained_upfront = False
        if return_fluxes and loopless == "fastSNP":
            add_loopless(model, method=loopless, reactions=cyclic_reactions)
            constrained_upfront = True

        for loopless_reactions, opt_rxn_ids in enumerate(reaction_ids_by_type):
            if len(opt_rxn_ids["minimum"]) == 0 and len(opt_rxn_ids["maximum"]) == 0:
                continue

            run_cycle_free_flux = bool(loopless_reactions)
            if loopless_reactions and loopless == "fastSNP":
                if not constrained_upfront:
                    add_loopless(
                        model,
                        method=loopless,
                        reactions=cyclic_reactions,
                    )
                run_cycle_free_flux = False

            for what in ("minimum", "maximum"):
                if len(opt_rxn_ids[what]) == 0:
                    continue

                cur_processes = min(processes, len(opt_rxn_ids[what]))
                if cur_processes > 1:
                    # We create and destroy a new pool here in order to set the
                    # objective direction for all reactions. This creates a
                    # slight overhead but seems the most clean.
                    chunk_size = len(opt_rxn_ids[what]) // cur_processes
                    with ProcessPool(
                        cur_processes,
                        initializer=_init_worker,
                        initargs=(
                            model,
                            run_cycle_free_flux,
                            what[:3],
                            return_fluxes,
                            time_limit,
                            accept_incumbent,
                        ),
                    ) as pool:
                        for rxn_id, value, fluxes, proven in pool.imap_unordered(
                            _fva_step, opt_rxn_ids[what], chunksize=chunk_size
                        ):
                            fva_result.at[rxn_id, what] = value
                            fva_result.at[rxn_id, what + "_proven"] = proven
                            if return_fluxes and fluxes is not None:
                                flux_rows[what][rxn_id] = fluxes
                else:
                    _init_worker(
                        model,
                        run_cycle_free_flux,
                        what[:3],
                        return_fluxes,
                        time_limit,
                        accept_incumbent,
                    )
                    for rxn_id, value, fluxes, proven in map(
                        _fva_step, opt_rxn_ids[what]
                    ):
                        fva_result.at[rxn_id, what] = value
                        fva_result.at[rxn_id, what + "_proven"] = proven
                        if return_fluxes and fluxes is not None:
                            flux_rows[what][rxn_id] = fluxes

    if time_limit is not None and processes <= 1:
        model.solver.configuration.timeout = orig_timeout

    if return_fluxes:
        return (
            fva_result[result_columns],
            {
                what: pd.DataFrame.from_dict(rows, orient="index")
                for what, rows in flux_rows.items()
            },
        )
    return fva_result[result_columns]


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
