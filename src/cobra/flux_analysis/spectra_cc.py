"""Provide an implementation of the SPECTRA consistency check."""

from logging import getLogger
from typing import TYPE_CHECKING, Dict, Optional, Set, Tuple

import numpy as np
from optlang.interface import OPTIMAL
from optlang.symbolics import Zero

from ..core import get_solution
from .helpers import normalize_cutoff


if TYPE_CHECKING:
    from cobra.core import Model, Reaction


logger = getLogger(__name__)

LARGE_VALUE = 1.0e6
#: Lower and upper limit of the uniformly sampled objective coefficients.
WEIGHT_RANGE = (1.0, 1.1)


def _reaction_sign(reaction: "Reaction") -> float:
    """Return the orientation a reaction is examined in.

    Reactions that can only carry negative flux are examined in their reverse
    orientation, which is equivalent to negating their column in the
    stoichiometric matrix. Every other reaction is examined as is.

    Parameters
    ----------
    reaction: cobra.Reaction
        The reaction to orient.

    Returns
    -------
    float
        -1.0 if the reaction can only carry negative flux, 1.0 otherwise.

    """
    return -1.0 if reaction.upper_bound <= 0.0 else 1.0


def _is_reversible(reaction: "Reaction") -> bool:
    """Return whether a reaction can carry flux in its reverse orientation.

    Reversibility is evaluated after orienting the reaction as described in
    :func:`_reaction_sign`, so a reaction restricted to negative flux counts
    as irreversible.

    Parameters
    ----------
    reaction: cobra.Reaction
        The reaction to inspect.

    Returns
    -------
    bool
        Whether the oriented reaction has a negative lower bound.

    """
    sign = _reaction_sign(reaction)
    return min(sign * reaction.lower_bound, sign * reaction.upper_bound) < 0.0


def _add_auxiliary_vars(
    model: "Model", signs: Dict[str, float], flux_threshold: float
) -> Set[str]:
    """Add the auxiliary variables and constraints for both LPs.

    Two families are added. The forward family bounds an auxiliary variable
    from above by ``flux_threshold`` and by the oriented flux, so that
    maximizing it drives the reaction towards positive flux. The reverse
    family does the mirror image for reversible reactions. Both are added in
    a disabled state; :func:`_find_flux_mode` enables the subset it needs.

    Parameters
    ----------
    model: cobra.Model
        The model to operate on.
    signs: dict of {str: float}
        The orientation of every reaction, keyed by reaction identifier.
    flux_threshold: float
        The magnitude the auxiliary variables, and hence the reaction fluxes,
        are driven to.

    Returns
    -------
    set of str
        The identifiers of the reversible reactions, for which a reverse
        auxiliary variable was added.

    """
    prob = model.problem
    reversible_ids = set()
    vars_and_cons = []

    for rxn in model.reactions:
        oriented_flux = signs[rxn.id] * rxn.flux_expression

        var = prob.Variable(
            f"spectra_fwd_aux_{rxn.id}", lb=-LARGE_VALUE, ub=flux_threshold
        )
        # Enabled by raising the lower bound back to 0, i.e. v_i >= z_i.
        const = prob.Constraint(
            oriented_flux - var, name=f"spectra_fwd_{rxn.id}", lb=-LARGE_VALUE
        )
        vars_and_cons.extend([var, const])

        if not _is_reversible(rxn):
            continue

        reversible_ids.add(rxn.id)
        var = prob.Variable(
            f"spectra_rev_aux_{rxn.id}", lb=-flux_threshold, ub=LARGE_VALUE
        )
        # Enabled by lowering the upper bound back to 0, i.e. v_i <= z_i.
        const = prob.Constraint(
            oriented_flux - var, name=f"spectra_rev_{rxn.id}", ub=LARGE_VALUE
        )
        vars_and_cons.extend([var, const])

    model.add_cons_vars(vars_and_cons)
    model.solver.update()
    return reversible_ids


def _find_flux_mode(
    model: "Model",
    rxn_ids: Set[str],
    zero_cutoff: float,
    rng: np.random.Generator,
    reverse: bool = False,
) -> Set[str]:
    """Perform one of the two LPs required for SPECTRA.

    The LP pushes as many of the given reactions as possible towards the flux
    threshold at once, in the forward or the reverse orientation. Randomized
    objective coefficients break ties between reactions, so that repeated
    calls explore different corners of the optimal face.

    Parameters
    ----------
    model: cobra.Model
        The model to operate on, already holding the auxiliary variables.
    rxn_ids: set of str
        The reactions to drive towards the flux threshold.
    zero_cutoff: float
        The cutoff below which a flux is considered zero.
    rng: numpy.random.Generator
        The random generator supplying the objective coefficients.
    reverse: bool, optional
        Whether to use the reverse family of auxiliary variables, driving the
        reactions towards negative flux instead (default False).

    Returns
    -------
    set of str
        The identifiers of the reactions carrying flux above the cutoff.

    """
    if not rxn_ids:
        return set()

    prefix = "spectra_rev" if reverse else "spectra_fwd"
    constraints = [model.constraints.get(f"{prefix}_{rid}") for rid in rxn_ids]

    for const in constraints:
        if reverse:
            const.ub = 0.0
        else:
            const.lb = 0.0

    try:
        weights = rng.uniform(*WEIGHT_RANGE, size=len(rxn_ids))
        obj_vars = [model.variables.get(f"{prefix}_aux_{rid}") for rid in rxn_ids]
        model.objective = model.problem.Objective(
            Zero, direction="min" if reverse else "max"
        )
        model.objective.set_linear_coefficients(dict(zip(obj_vars, weights)))

        model.slim_optimize()
        status = model.solver.status
        if status == OPTIMAL:
            fluxes = get_solution(model).fluxes
            return set(fluxes[fluxes.abs() >= zero_cutoff].index)

        logger.warning(
            "The %s LP terminated with status '%s'; treating its %d reactions "
            "as carrying no flux.",
            "reverse" if reverse else "forward",
            status,
            len(rxn_ids),
        )
        return set()
    finally:
        # Leave the constraints slack again, so that the next LP starts from a
        # clean problem even if this one failed.
        for const in constraints:
            if reverse:
                const.ub = LARGE_VALUE
            else:
                const.lb = -LARGE_VALUE


def _find_consistent_reaction_ids(
    model: "Model",
    tol: float = 1e-4,
    zero_cutoff: Optional[float] = None,
    seed: Optional[int] = None,
) -> Tuple[Set[str], int]:
    """Find the flux consistent reactions of a model using SPECTRA [1]_.

    Parameters
    ----------
    model: cobra.Model
        The model to operate on.
    tol: float, optional
        The flux threshold to consider, i.e. the magnitude the LPs drive the
        reaction fluxes to (default 1e-4).
    zero_cutoff: float, optional
        The cutoff below which a flux is considered zero (default
        `model.tolerance`).
    seed: int, optional
        A seed for the random objective coefficients, making the result
        reproducible (default None).

    Returns
    -------
    tuple of (set of str, int)
        The identifiers of the consistent reactions and the number of LPs
        that were solved to find them.

    References
    ----------
    .. [1] S, P. K., Sridhar, S., Alsmadi, N., Mahadevan, R., & Bhatt, N. P.
           (2026). Generalist method to reconstruct metabolic networks from
           multi-omics data at large-scale. bioRxiv.
           https://doi.org/10.64898/2026.04.02.716249

    """
    zero_cutoff = normalize_cutoff(model, zero_cutoff)
    rng = np.random.default_rng(seed)

    signs = {rxn.id: _reaction_sign(rxn) for rxn in model.reactions}
    # Reactions that have not been shown to carry flux yet.
    rxns_to_check = set(signs)
    n_lps = 0

    with model:
        reversible_ids = _add_auxiliary_vars(model, signs, tol)

        previous_count = None
        while len(rxns_to_check) != previous_count:
            previous_count = len(rxns_to_check)
            n_lps += 2

            rxns_to_check -= _find_flux_mode(
                model, rxns_to_check, zero_cutoff, rng, reverse=False
            )
            rxns_to_check -= _find_flux_mode(
                model,
                rxns_to_check & reversible_ids,
                zero_cutoff,
                rng,
                reverse=True,
            )
            logger.debug(
                "LPs solved: %d - reactions left to check: %d",
                n_lps,
                len(rxns_to_check),
            )

    consistent_ids = set(signs) - rxns_to_check
    logger.info(
        "Final - consistent reactions: %d - inconsistent reactions: %d "
        "[LPs=%d, tol=%.2g, cutoff=%.2g]",
        len(consistent_ids),
        len(rxns_to_check),
        n_lps,
        tol,
        zero_cutoff,
    )
    return consistent_ids, n_lps


def spectra_cc(
    model: "Model",
    tol: float = 1e-4,
    zero_cutoff: Optional[float] = None,
    seed: Optional[int] = None,
) -> "Model":
    r"""
    Check consistency of a metabolic network using SPECTRA [1]_.

    SPECTRA's consistency check is a pure LP method for removing the blocked
    reactions of a metabolic network. Like FASTCC it drives as many reactions
    as possible towards a flux threshold at once, but it handles reversible
    reactions without flipping them: the auxiliary variable of the forward LP
    is unbounded from below, so a reaction that can only carry negative flux
    simply does not contribute to that LP, and is picked up by the reverse LP
    instead. The two LPs alternate until the set of unexplained reactions
    stops shrinking. For more details, please check [1]_.

    Parameters
    ----------
    model: cobra.Model
        The model to operate on.
    tol: float, optional
        The flux threshold to consider, i.e. the magnitude the LPs drive the
        reaction fluxes to (default 1e-4).
    zero_cutoff: float, optional
        The cutoff below which a flux is considered zero (default
        `model.tolerance`).
    seed: int, optional
        A seed for the random objective coefficients, making the result
        reproducible (default None).

    Returns
    -------
    cobra.Model
        The consistent model.

    Notes
    -----
    The pair of LPs used for SPECTRA is like so:
    maximize: \sum_{i \in J} w_i z_i
    s.t.    : z_i \in [-\infty, \varepsilon] \forall i \in J
              v_i \ge z_i \forall i \in J
              Sv = 0, v \in B

    minimize: \sum_{i \in J^{rev}} w_i z_i
    s.t.    : z_i \in [-\varepsilon, \infty] \forall i \in J^{rev}
              v_i \le z_i \forall i \in J^{rev}
              Sv = 0, v \in B

    where :math:`J` are the reactions not yet shown to carry flux,
    :math:`J^{rev}` the reversible ones among them, and :math:`w_i` are
    coefficients drawn uniformly from [1, 1.1].

    A reaction is taken to be consistent as soon as it carries a flux above
    `zero_cutoff` in any of these LPs, which is the same criterion
    :func:`~cobra.flux_analysis.fastcc.fastcc` applies, so both functions
    agree on which reactions are blocked. The reference MATLAB implementation
    instead requires a reaction to reach the flux threshold itself; pass
    ``zero_cutoff=0.99 * tol`` to reproduce that, at the cost of also
    discarding reactions that can carry some flux but never as much as `tol`.

    Because a reaction that resists the forward LP is handled by a dedicated
    reverse LP instead of being retried, :func:`spectra_cc` never needs
    :func:`~cobra.flux_analysis.fastcc.fastcc`'s per-reaction fallback (its
    ``singletons`` phase). On models with few reversible reactions the two
    take a similar number of LP solves; on larger, more reversibility-rich
    models :func:`spectra_cc` tends to need markedly fewer, while returning
    the same consistent set. See ``benchmarks/spectra_cc_vs_fastcc.ipynb``
    for a runtime comparison on bundled test models.

    References
    ----------
    .. [1] S, P. K., Sridhar, S., Alsmadi, N., Mahadevan, R., & Bhatt, N. P.
           (2026). Generalist method to reconstruct metabolic networks from
           multi-omics data at large-scale. bioRxiv.
           https://doi.org/10.64898/2026.04.02.716249

    """
    consistent_ids, _ = _find_consistent_reaction_ids(model, tol, zero_cutoff, seed)

    consistent_model = model.copy()
    consistent_model.remove_reactions(
        {rxn.id for rxn in model.reactions} - consistent_ids, remove_orphans=True
    )

    return consistent_model
