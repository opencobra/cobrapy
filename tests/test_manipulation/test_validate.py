"""Test functionalities of cobra component validation functions."""

import pytest

from cobra.core import Metabolite, Model, Reaction
from cobra.manipulation import check_mass_balance, check_metabolite_compartment_formula


def test_validate_mass_balance(model: Model) -> None:
    """Test reaction mass balance validation."""
    assert len(check_mass_balance(model)) == 0
    # if we remove the SBO term which marks the reaction as
    # mass balanced, then the reaction should be detected as
    # no longer mass balanced
    EX_rxn = model.reactions.query(lambda r: r.boundary)[0]
    EX_rxn.annotation.pop("sbo")
    balance = check_mass_balance(model)
    assert len(balance) == 1
    assert EX_rxn in balance
    m1 = Metabolite("m1", formula="()")
    r1 = Reaction("r1")
    r1.add_metabolites({m1: 1})
    with pytest.raises(ValueError), pytest.warns(UserWarning):
        r1.check_mass_balance()


def test_validate_formula_compartment(model: Model) -> None:
    """Test metabolite formulae validation."""
    model.metabolites[1].formula = "(a*.bcde)"
    errors = check_metabolite_compartment_formula(model)
    assert len(errors) == 1
