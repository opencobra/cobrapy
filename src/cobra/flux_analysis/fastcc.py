"""Provide an implementation of FASTCC."""

from logging import getLogger
from typing import TYPE_CHECKING, List, Optional, Set

from optlang.symbolics import Zero

from .helpers import normalize_cutoff


if TYPE_CHECKING:
    from cobra.core import Model, Reaction

logger = getLogger(__name__)
LARGE_VALUE = 1.0e6


def _add_lp7_vars(
    model: "Model", rxns: List["Reaction"], flux_threshold: float
) -> None:
    """Add the variables and constraints for the LP.

    Parameters
    ----------
    model: cobra.Model
        The model to operate on.
    rxns: list of cobra.Reaction
        The reactions to use for LP.
    flux_threshold: float
        The upper threshold an auxiliary variable can have.

    """
    prob = model.problem
    vars_and_cons = []

    for rxn in rxns:
        var = prob.Variable(f"auxiliary_{rxn.id}", lb=0.0, ub=flux_threshold)
        const = prob.Constraint(
            rxn.flux_expression - var,
            name="aux_constraint_{}".format(rxn.id),
            lb=-LARGE_VALUE,
        )
        vars_and_cons.extend([var, const])

    model.add_cons_vars(vars_and_cons)
    model.solver.update()


def _find_sparse_mode(
    model: "Model", rxn_ids: Set[str], zero_cutoff: float
) -> List["Reaction"]:
    """Perform the LP required for FASTCC.

    Parameters
    ----------
    model: cobra.Model
        The model to perform FASTCC on.
    rxns: list of cobra.Reaction
        The reactions to use for LP.
    zero_cutoff: float
        The cutoff below which flux is considered zero.

    Returns
    -------
    list of cobra.Reaction
        The list of reactions to consider as consistent.

    """
    if not rxn_ids:
        return set()

    # Enable constraints for the reactions
    for rid in rxn_ids:
        model.constraints.get(f"aux_constraint_{rid}").lb = 0.0

    obj_vars = [model.variables.get(f"auxiliary_{rid}") for rid in rxn_ids]
    model.objective = Zero
    model.objective.set_linear_coefficients({v: 1.0 for v in obj_vars})

    sol = model.optimize(objective_sense="max")

    # Disable constraints for the reactions
    for rid in rxn_ids:
        model.constraints.get("aux_constraint_{}".format(rid)).lb = -LARGE_VALUE

    return set(sol.fluxes[sol.fluxes.abs() > zero_cutoff].index)


def _flip_coefficients(model: "Model", rxn_ids: Set[str]) -> None:
    """Flip the coefficients for optimizing in reverse direction.

    Parameters
    ----------
    model: cobra.Model
        The model to operate on.
    rxns: list of cobra.Reaction
        The list of reactions whose coefficients will be flipped.

    """
    if not rxn_ids:
        return
    # flip reactions
    for rxn in rxn_ids:
        const = model.constraints.get(f"aux_constraint_{rxn}")
        var = model.variables.get(f"auxiliary_{rxn}")
        coefs = const.get_linear_coefficients(const.variables)
        const.set_linear_coefficients({k: -v for k, v in coefs.items() if k is not var})
        model.solver.update()


def _any_set(s):
    for x in s:
        return set([x])


def fastcc(
    model: "Model", flux_threshold: float = 1.0, zero_cutoff: Optional[float] = None
) -> "Model":
    r"""
    Check consistency of a metabolic network using FASTCC [1]_.

    FASTCC (Fast Consistency Check) is an algorithm for rapid and
    efficient consistency check in metabolic networks. FASTCC is
    a pure LP implementation and is low on computation resource
    demand. FASTCC also circumvents the problem associated with
    reversible reactions for the purpose. Given a global model,
    it will generate a consistent global model i.e., remove
    blocked reactions. For more details on FASTCC, please
    check [1]_.

    Parameters
    ----------
    model: cobra.Model
        The model to operate on.
    flux_threshold: float, optional
        The flux threshold to consider (default 1.0).
    zero_cutoff: float, optional
        The cutoff to consider for zero flux (default model.tolerance).

    Returns
    -------
    cobra.Model
        The consistent model.

    Notes
    -----
    The LP used for FASTCC is like so:
    maximize: \sum_{i \in J} z_i
    s.t.    : z_i \in [0, \varepsilon] \forall i \in J, z_i \in \mathbb{R}_+
              v_i \ge z_i \forall i \in J
              Sv = 0 v \in B

    References
    ----------
    .. [1] Vlassis N, Pacheco MP, Sauter T (2014)
           Fast Reconstruction of Compact Context-Specific Metabolic Network
           Models.
           PLoS Comput Biol 10(1): e1003424. doi:10.1371/journal.pcbi.1003424

    """
    zero_cutoff = normalize_cutoff(model, zero_cutoff)

    all_rxns = {rxn.id for rxn in model.reactions}
    irreversible_rxns = {rxn.id for rxn in model.reactions if not rxn.reversibility}
    rxns_to_check = irreversible_rxns
    flipped = False
    singletons = False

    with model:
        _add_lp7_vars(model, model.reactions, flux_threshold)

        rxns_to_keep = _find_sparse_mode(model, rxns_to_check, zero_cutoff)
        rxns_to_check = all_rxns.difference(rxns_to_keep)
        logger.info(
            "Initial step found %d consistent reactions. "
            "Starting the consistency loop for the remaining %d reactions.",
            len(rxns_to_keep),
            len(rxns_to_check),
        )

        while rxns_to_check:
            logger.debug(
                "reactions to check: %d - consistent reactions:"
                " %d - flipped: %d - singletons: %d",
                len(rxns_to_check),
                len(rxns_to_keep),
                flipped,
                singletons,
            )
            check = _any_set(rxns_to_check) if singletons else rxns_to_check
            new_rxns = _find_sparse_mode(model, check, zero_cutoff)
            rxns_to_keep.update(new_rxns)

            if rxns_to_check.intersection(rxns_to_keep):
                rxns_to_check = rxns_to_check.difference(rxns_to_keep)
                flipped = False
            else:
                check_irr = check.difference(irreversible_rxns)
                if flipped or not check_irr:
                    if singletons:
                        logger.debug("%s is inconsistent", check)
                        rxns_to_check = rxns_to_check.difference(check)
                    flipped = False
                    singletons = True
                else:
                    flipped = True
                    check = _any_set(rxns_to_check) if singletons else rxns_to_check
                    _flip_coefficients(model, check_irr)
        logger.info(
            "Final - consistent reactions: %d"
            " - inconsistent reactions: %d [eps=%.2g, tol=%.2g]",
            len(rxns_to_keep),
            len(all_rxns) - len(rxns_to_keep),
            flux_threshold,
            zero_cutoff,
        )

    consistent_model = model.copy()
    consistent_model.remove_reactions(
        all_rxns.difference(rxns_to_keep), remove_orphans=True
    )

    return consistent_model
