"""Provide an implementation of the SPECTRA consistency check."""

from copy import deepcopy
from logging import getLogger
from typing import TYPE_CHECKING, Dict, Iterable, Optional, Set, Tuple

import numpy as np
from optlang.interface import OPTIMAL
from optlang.symbolics import Zero

from ..core import Group, Model
from ..core.dictlist import DictList
from .helpers import normalize_cutoff


if TYPE_CHECKING:
    from cobra.core import Reaction


logger = getLogger(__name__)

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


def _find_flux_mode(
    model: "Model",
    rxn_ids: Iterable[str],
    signs: Dict[str, float],
    flux_threshold: float,
    zero_cutoff: float,
    rng: np.random.Generator,
    reverse: bool = False,
) -> Set[str]:
    """Perform one of the two LPs required for SPECTRA.

    An auxiliary variable is tied to each given reaction's oriented flux by a
    one-sided constraint, and a randomly weighted sum of those auxiliary
    variables is pushed as far as the flux bounds allow. The randomized
    weights break ties between reactions, so that repeated calls explore
    different corners of the optimal face.

    The auxiliary variables are free on the side they are pushed away from,
    rather than held by a large finite bound. A finite bound would quietly
    become binding on a model whose flux bounds exceed it.

    Parameters
    ----------
    model: cobra.Model
        The model to operate on.
    rxn_ids: iterable of str
        The reactions to drive towards the flux threshold.
    signs: dict of {str: float}
        The orientation of every reaction, keyed by reaction identifier.
    flux_threshold: float
        The magnitude the auxiliary variables, and hence the reaction fluxes,
        are driven to.
    zero_cutoff: float
        The cutoff below which a flux is considered zero.
    rng: numpy.random.Generator
        The random generator supplying the objective coefficients.
    reverse: bool, optional
        Whether to drive the reactions towards negative oriented flux instead
        (default False).

    Returns
    -------
    set of str
        The identifiers of the reactions carrying flux above the cutoff.

    """
    # Sorted, not merely listed: the objective coefficients are drawn in one
    # block and zipped against the variables positionally, so the order
    # decides which reaction gets which weight. Callers pass a set, whose
    # iteration order varies with PYTHONHASHSEED, which would make a seeded
    # run irreproducible between processes.
    rxn_ids = sorted(rxn_ids)
    if not rxn_ids:
        return set()

    prob = model.problem
    aux_vars = []
    cons_and_vars = []
    coefficients = []

    for rxn_id in rxn_ids:
        rxn = model.reactions.get_by_id(rxn_id)
        sign = signs[rxn_id]

        if reverse:
            # z_i >= -eps and v_i <= z_i, minimized: drives v_i down to -eps.
            var = prob.Variable(f"spectra_aux_{rxn_id}", lb=-flux_threshold)
            const = prob.Constraint(Zero, name=f"spectra_cons_{rxn_id}", ub=0.0)
        else:
            # z_i <= eps and v_i >= z_i, maximized: drives v_i up to eps.
            var = prob.Variable(f"spectra_aux_{rxn_id}", ub=flux_threshold)
            const = prob.Constraint(Zero, name=f"spectra_cons_{rxn_id}", lb=0.0)

        aux_vars.append(var)
        cons_and_vars.extend([var, const])
        # The constraint is sign_i * v_i - z_i, and v_i is itself the
        # difference of the forward and reverse variables. Declaring it empty
        # and filling in the three coefficients avoids handing optlang an
        # expression to build symbolically, which sympy otherwise dominates
        # the runtime of.
        coefficients.append(
            (
                const,
                {
                    rxn.forward_variable: sign,
                    rxn.reverse_variable: -sign,
                    var: -1.0,
                },
            )
        )

    with model:
        model.add_cons_vars(cons_and_vars)
        # Coefficients can only be set once the constraints belong to a
        # problem.
        model.solver.update()
        for const, terms in coefficients:
            const.set_linear_coefficients(terms)
        # Reversed in place so that leaving the context deletes from the tail
        # of optlang's containers: removing from the front re-indexes
        # everything after it, which is quadratic in the number of auxiliary
        # variables. The context holds this very list, so reversing it here
        # reverses the teardown.
        cons_and_vars.reverse()

        model.objective = prob.Objective(Zero, direction="min" if reverse else "max")
        weights = rng.uniform(*WEIGHT_RANGE, size=len(aux_vars))
        model.objective.set_linear_coefficients(dict(zip(aux_vars, weights)))

        model.slim_optimize()
        status = model.solver.status
        if status != OPTIMAL:
            logger.warning(
                "The %s LP terminated with status '%s'; treating its %d "
                "reactions as carrying no flux.",
                "reverse" if reverse else "forward",
                status,
                len(rxn_ids),
            )
            return set()

        return {rxn.id for rxn in model.reactions if abs(rxn.flux) >= zero_cutoff}


def _subset_model(model: "Model", keep_ids: Set[str]) -> "Model":
    """Return a new model holding only the given reactions.

    The model is *built* from the reactions that are kept rather than copied
    whole and then cut down. Deleting is quadratic: each removal makes
    ``optlang.Container._reindex`` rebuild the whole name-to-index map, and
    the solver is updated once per reaction. Building instead is linear in
    what is kept, and populates the solver in a single batch.

    Parameters
    ----------
    model: cobra.Model
        The model to take the subset of. It is left unchanged.
    keep_ids: set of str
        The identifiers of the reactions to keep.

    Returns
    -------
    cobra.Model
        A new model carrying the kept reactions and everything attached to
        them, in the order the source had them.

    """
    keep = set(keep_ids)
    reactions = [rxn for rxn in model.reactions if rxn.id in keep]
    metabolites = {met for rxn in reactions for met in rxn._metabolites}
    genes = {gene for rxn in reactions for gene in rxn._genes}

    # ``deepcopy`` carries every attribute, so nothing has to be enumerated
    # and nothing is silently dropped. It needs ``_model`` detached first:
    # ``Reaction.__getstate__`` does not null it, unlike
    # ``Object.__getstate__``, so a plain deepcopy would drag the whole source
    # model, and its solver, along behind each reaction.
    detached = [(obj, obj._model) for obj in (*reactions, *metabolites, *genes)]
    for obj, _ in detached:
        obj._model = None
    try:
        copies = deepcopy(reactions)
    finally:
        for obj, owner in detached:
            obj._model = owner

    subset = Model(model.id)
    subset.name = model.name
    subset.compartments = dict(model.compartments)
    subset.annotation = deepcopy(model.annotation)
    subset.notes = deepcopy(model.notes)
    subset.add_reactions(copies)

    # ``Reaction.__getstate__`` serialises the GPR to a string, so the rule is
    # re-parsed when the reactions are added and the genes come back as fresh
    # objects carrying nothing but an identifier. Their metadata has to be put
    # back from the source, or a model whose genes are annotated loses every
    # annotation it had.
    for gene in subset.genes:
        original = model.genes.get_by_id(gene.id)
        gene.name = original.name
        gene.annotation = deepcopy(original.annotation)
        gene.notes = deepcopy(original.notes)

    # ``add_reactions`` interns metabolites and genes in the order it first
    # meets them, which is not the order they had in the source.
    for attribute in ("metabolites", "genes"):
        position = {obj.id: i for i, obj in enumerate(getattr(model, attribute))}
        current = list(getattr(subset, attribute))
        setattr(
            subset,
            attribute,
            DictList(sorted(current, key=lambda obj: position[obj.id])),
        )

    # Groups survive the loss of all their members, as they would have done
    # when the reactions were deleted out from under them.
    groups = []
    for group in model.groups:
        members = []
        for member in group.members:
            for container in (subset.reactions, subset.metabolites, subset.genes):
                if member.id in container:
                    members.append(container.get_by_id(member.id))
                    break
        new_group = Group(group.id, name=group.name, members=members, kind=group.kind)
        new_group.annotation = deepcopy(group.annotation)
        new_group.notes = deepcopy(group.notes)
        groups.append(new_group)
    if groups:
        subset.add_groups(groups)

    # Nothing above reaches model-level state, so it is carried over by hand.
    subset.solver = model.solver.interface
    subset.tolerance = model.tolerance
    coefficients = {
        subset.reactions.get_by_id(rxn.id): rxn.objective_coefficient
        for rxn in reactions
        if rxn.objective_coefficient
    }
    if coefficients:
        subset.objective = coefficients
    # Unconditional: with every objective reaction gone the coefficients are
    # empty but the direction still has to follow the original.
    subset.objective.direction = model.objective.direction
    return subset


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
    reversible_ids = {rxn.id for rxn in model.reactions if _is_reversible(rxn)}
    # Reactions that have not been shown to carry flux yet.
    rxns_to_check = set(signs)
    n_lps = 0

    previous_count = None
    while len(rxns_to_check) != previous_count:
        previous_count = len(rxns_to_check)
        n_lps += 2

        rxns_to_check -= _find_flux_mode(
            model, rxns_to_check, signs, tol, zero_cutoff, rng, reverse=False
        )
        rxns_to_check -= _find_flux_mode(
            model,
            rxns_to_check & reversible_ids,
            signs,
            tol,
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

    References
    ----------
    .. [1] S, P. K., Sridhar, S., Alsmadi, N., Mahadevan, R., & Bhatt, N. P.
           (2026). Generalist method to reconstruct metabolic networks from
           multi-omics data at large-scale. bioRxiv.
           https://doi.org/10.64898/2026.04.02.716249

    """
    consistent_ids, _ = _find_consistent_reaction_ids(model, tol, zero_cutoff, seed)

    return _subset_model(model, consistent_ids)
