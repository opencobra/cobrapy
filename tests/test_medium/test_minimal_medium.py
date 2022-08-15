"""Test functionalities of minimal medium creation and analysis."""

from typing import Callable

import pandas as pd
import pytest

from cobra.core import Model, Reaction
from cobra.medium import minimal_medium


def test_minimal_medium_linear(model: Model) -> None:
    """Test linear minimal medium."""
    med = minimal_medium(model, 0.8)
    assert len(med) <= 4
    assert all(med > 1e-6)


def test_minimal_medium_linear_benchmark(model: Model, benchmark: Callable) -> None:
    """Benchmark linear minimal medium."""
    benchmark(minimal_medium, model, 0.8)


def test_minimal_medium_mip(model: Model) -> None:
    """Test mixed-integer minimal medium."""
    med = minimal_medium(model, 0.8, minimize_components=True)
    assert len(med) <= 4
    assert all(med > 1e-6)
    # Anaerobic growth
    med = minimal_medium(model, 0.1, minimize_components=True)
    assert len(med) <= 3
    assert all(med > 1e-6)


def test_minimal_medium_mip_benchmark(model: Model, benchmark: Callable) -> None:
    """Benchmark mixed-integer minimal medium."""
    benchmark(minimal_medium, model, 0.8, True)


def test_minimal_medium_alternative_mip(model: Model) -> None:
    """Test alternative mixed-integer minimal medium."""
    med = minimal_medium(model, 0.8, minimize_components=5, open_exchanges=True)
    assert isinstance(med, pd.DataFrame)
    assert med.shape[0] >= 5
    assert med.shape[1] == 5
    assert all((med > 0).sum() == 3)
    assert all(med.sum(axis=1) > 1e-6)


def test_minimal_medium_exports(model: Model) -> None:
    """Test exports of a minimal medium."""
    med = minimal_medium(model, 0.8, exports=True, minimize_components=True)
    assert len(med) > 4
    assert any(med < -1e-6)


def test_minimal_medium_open_exchanges(model: Model) -> None:
    """Test open exchanges of a minimal medium."""
    model.reactions.EX_glc__D_e.bounds = 0, 0
    med = minimal_medium(model, 0.8)
    assert med is None
    med = minimal_medium(model, 0.8, minimize_components=True)
    assert med is None

    med = minimal_medium(model, 0.8, open_exchanges=True)
    assert len(med) >= 3
    med = minimal_medium(model, 0.8, open_exchanges=100)
    assert len(med) >= 3


def test_model_medium(model: Model) -> None:
    """Test proper functioning of model medium manipulations."""
    # Add a dummy 'malformed' import reaction
    bad_import = Reaction("bad_import")
    bad_import.add_metabolites({model.metabolites.pyr_c: 1})
    bad_import.bounds = (0, 42)
    model.add_reactions([bad_import])

    # Test basic setting and getting methods
    medium = model.medium
    model.medium = medium
    assert model.medium == medium

    # Test context management
    with model:
        # Ensure the bounds are correct beforehand
        assert model.reactions.EX_glc__D_e.lower_bound == -10
        assert model.reactions.bad_import.upper_bound == 42
        assert model.reactions.EX_co2_e.lower_bound == -1000

        # Make changes to the media
        new_medium = model.medium
        new_medium["EX_glc__D_e"] = 20
        new_medium["bad_import"] = 24
        del new_medium["EX_co2_e"]

        # Change the medium, make sure changes work
        model.medium = new_medium
        assert model.reactions.EX_glc__D_e.lower_bound == -20
        assert model.reactions.bad_import.upper_bound == 24
        assert model.reactions.EX_co2_e.lower_bound == 0

    # Make sure changes revert after the contex
    assert model.reactions.EX_glc__D_e.lower_bound == -10
    assert model.reactions.bad_import.upper_bound == 42
    assert model.reactions.EX_co2_e.lower_bound == -1000

    new_medium["bogus_rxn"] = 0
    with pytest.raises(KeyError):
        model.medium = new_medium


def test_medium_does_not_affect_reactant_exports(model: Model) -> None:
    """Test that the medium setter does not overwrite exports defined as reactants."""
    med = model.medium
    # Set a fixed export rate
    model.reactions.EX_ac_e.lower_bound = 0.1
    model.medium = med
    assert model.reactions.EX_ac_e.lower_bound == 0.1

    # should be overwritten if actually in the medium
    med["EX_ac_e"] = 1
    model.medium = med
    assert model.reactions.EX_ac_e.lower_bound == -1


def test_medium_does_not_affect_product_exports(model: Model) -> None:
    """Test that the medium setter does not overwrite exports defined as products."""
    med = model.medium
    # Reverse reaction definition
    model.reactions.EX_ac_e *= -1
    # Set a fixed export rate
    model.reactions.EX_ac_e.bounds = -1, -0.1
    model.medium = med
    assert model.reactions.EX_ac_e.bounds == (-1, -0.1)

    # should be overwritten if actually in the medium
    med["EX_ac_e"] = 1
    model.medium = med
    assert model.reactions.EX_ac_e.upper_bound == 1
