"""Test functionalities of Phenotype Phase Plane Analysis."""

import numpy as np
import pytest

from cobra import Model
from cobra.flux_analysis import production_envelope


def test_envelope_one(model: Model) -> None:
    """Test flux of production envelope."""
    df = production_envelope(model, ["EX_o2_e"])
    assert np.isclose(df["flux_maximum"].sum(), 9.342, atol=1e-3)


def test_envelope_multi_reaction_objective(model: Model) -> None:
    """Test production of multiple objectives."""
    obj = {model.reactions.EX_ac_e: 1, model.reactions.EX_co2_e: 1}
    with pytest.raises(ValueError):
        production_envelope(model, "EX_o2_e", obj)


@pytest.mark.parametrize(
    "variables, num",
    [
        (["EX_glc__D_e"], 30),
        (["EX_glc__D_e", "EX_o2_e"], 20),
        (["EX_glc__D_e", "EX_o2_e", "EX_ac_e"], 10),
    ],
)
def test_multi_variable_envelope(model: Model, variables: str, num: int) -> None:
    """Test production of envelope (multiple variable)."""
    df = production_envelope(model, variables, points=num)
    assert len(df) == num ** len(variables)


def test_envelope_two(model: Model) -> None:
    """Test production of envelope."""
    df = production_envelope(model, ["EX_glc__D_e", "EX_o2_e"], objective="EX_ac_e")
    assert np.isclose(df["flux_maximum"].sum(), 1737.466, atol=1e-3)
    assert np.isclose(df["carbon_yield_maximum"].sum(), 83.579, atol=1e-3)
    assert np.isclose(df["mass_yield_maximum"].sum(), 82.176, atol=1e-3)
