"""Provide an implementation of FASTCC."""

from typing import TYPE_CHECKING, List, Optional

from optlang.symbolics import Zero

from .helpers import normalize_cutoff


if TYPE_CHECKING:
    from cobra.core import Model, Reaction


def _find_sparse_mode(
    model: "Model", rxns: List["Reaction"], flux_threshold: float, zero_cutoff: float
) -> List["Reaction"]:
    """Perform the LP required for FASTCC.

    Parameters
    ----------
    model: cobra.Model
        The model to perform FASTCC on.
    rxns: list of cobra.Reaction
        The reactions to use for LP.
    flux_threshold: float
        The upper threshold an auxiliary variable can have.
    zero_cutoff: float
        The cutoff below which flux is considered zero.

    Returns
    -------
    list of cobra.Reaction
        The list of reactions to consider as consistent.

    """
    if rxns:
        obj_vars = []
        vars_and_cons = []
        prob = model.problem

        for rxn in rxns:
            var = prob.Variable(
                "auxiliary_{}".format(rxn.id), lb=0.0, ub=flux_threshold
            )
            const = prob.Constraint(
                rxn.forward_variable + rxn.reverse_variable - var,
                name="constraint_{}".format(rxn.id),
                lb=0.0,
            )
            vars_and_cons.extend([var, const])
            obj_vars.append(var)

        model.add_cons_vars(vars_and_cons)
        model.objective = prob.Objective(Zero, sloppy=True)
        model.objective.set_linear_coefficients({v: 1.0 for v in obj_vars})

        model.optimize(objective_sense="max")
        result = [rxn for rxn in model.reactions if abs(rxn.flux) > zero_cutoff]
    else:
        result = []

    return result


def _flip_coefficients(model: "Model", rxns: List["Reaction"]) -> None:
    """Flip the coefficients for optimizing in reverse direction.

    Parameters
    ----------
    model: cobra.Model
        The model to operate on.
    rxns: list of cobra.Reaction
        The list of reactions whose coefficients will be flipped.

    """
    # flip reactions
    for rxn in rxns:
        const = model.constraints.get("constraint_{}".format(rxn.id))
        var = model.variables.get("auxiliary_{}".format(rxn.id))
        coefs = const.get_linear_coefficients(const.variables)
        const.set_linear_coefficients({k: -v for k, v in coefs.items() if k is not var})

    # flip objective
    objective = model.objective
    objective_coefs = objective.get_linear_coefficients(objective.variables)
    objective.set_linear_coefficients({k: -v for k, v in objective_coefs.items()})


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

    irreversible_rxns = [rxn for rxn in model.reactions if not rxn.reversibility]
    rxns_to_check = irreversible_rxns

    with model:
        rxns_to_keep = _find_sparse_mode(
            model, rxns_to_check, flux_threshold, zero_cutoff
        )

    rxns_to_check = list(set(model.reactions).difference(rxns_to_keep))

    while rxns_to_check:
        with model:
            new_rxns = _find_sparse_mode(
                model, rxns_to_check, flux_threshold, zero_cutoff
            )
            rxns_to_keep.extend(new_rxns)

            # this condition will be valid for all but the last iteration
            if list(set(rxns_to_check).intersection(rxns_to_keep)):
                rxns_to_check = list(set(rxns_to_check).difference(rxns_to_keep))

            else:
                rxns_to_flip = list(set(rxns_to_check).difference(irreversible_rxns))
                _flip_coefficients(model, rxns_to_flip)
                sol = model.optimize(min)
                to_add_rxns = sol.fluxes.index[sol.fluxes.abs() > zero_cutoff].tolist()
                rxns_to_keep.extend(
                    [model.reactions.get_by_id(rxn) for rxn in to_add_rxns]
                )
                # since this is the last iteration, it needs to break or else
                # it will run forever since rxns_to_check won't be empty
                break

    consistent_rxns = set(rxns_to_keep)
    # need the ids since Reaction objects are created fresh with model.copy()
    rxns_to_remove = [
        rxn.id for rxn in set(model.reactions).difference(consistent_rxns)
    ]

    consistent_model = model.copy()
    consistent_model.remove_reactions(rxns_to_remove, remove_orphans=True)

    return consistent_model
