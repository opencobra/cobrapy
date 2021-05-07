"""Provide functions for model component validations."""

from typing import TYPE_CHECKING, Dict, List


if TYPE_CHECKING:
    from cobra import Metabolite, Model, Reaction


# Set of mass unbalanced SBO terms
_NOT_MASS_BALANCED_TERMS = {
    "SBO:0000627",  # EXCHANGE
    "SBO:0000628",  # DEMAND
    "SBO:0000629",  # BIOMASS
    "SBO:0000631",  # PSEUDOREACTION
    "SBO:0000632",  # SINK
}


def check_mass_balance(model: "Model") -> Dict["Reaction", Dict["Metabolite", float]]:
    """Check mass balance for reactions of `model` and return unbalanced ones.

    Parameters
    ----------
    model: cobra.Model
        The model to perform check on.

    Returns
    -------
    dict of {cobra.Reaction: dict of {cobra.Metabolite: float}}
        Returns an empty dict if all components are balanced.

    """
    unbalanced = {}
    for reaction in model.reactions:
        if reaction.annotation.get("sbo") not in _NOT_MASS_BALANCED_TERMS:
            balance = reaction.check_mass_balance()
            if balance:
                unbalanced[reaction] = balance
    return unbalanced


def check_metabolite_compartment_formula(model: "Model") -> List[str]:
    """Check metabolite formulae of `model`.

    Parameters
    ----------
    model: cobra.Model
        The model to perform check on.

    Returns
    -------
    list of str
        Returns an empty list if no errors are found.

    """
    errors = []
    for met in model.metabolites:
        if met.formula is not None and len(met.formula) > 0:
            if not met.formula.isalnum():
                errors.append(
                    f"Metabolite '{met.id}' formula '{met.formula}' not alphanumeric"
                )
    return errors
