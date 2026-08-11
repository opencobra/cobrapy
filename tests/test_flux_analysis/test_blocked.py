"""Test functionalities for blocked reaction analysis."""

from typing import List

import pytest

from cobra import Model
from cobra.flux_analysis import find_blocked_reactions


def test_find_blocked_reactions_solver_none(model: Model) -> None:
    """Test find_blocked_reactions() [no specific solver]."""
    result = find_blocked_reactions(model, model.reactions[40:46])
    assert result == ["FRUpts2"]


def test_find_blocked_reactions(model: Model, all_solvers: List[str]) -> None:
    """Test find_blocked_reactions()."""
    model.solver = all_solvers
    result = find_blocked_reactions(model, model.reactions[40:46])
    assert result == ["FRUpts2"]

    result = find_blocked_reactions(model, model.reactions[42:48])
    assert set(result) == {"FUMt2_2", "FRUpts2"}

    result = find_blocked_reactions(model, model.reactions[30:50], open_exchanges=True)
    assert result == []


def test_find_blocked_reactions_cycle_free_flux_unsupported(model: Model) -> None:
    """Test find_blocked_reactions() rejects approximate loopless FVA."""
    with pytest.raises(ValueError, match="cycleFreeFlux"):
        find_blocked_reactions(model, loopless="cycleFreeFlux")


@pytest.mark.parametrize("loopless", ["fastSNP", "potentials"])
def test_find_blocked_reactions_with_loopless(
    model: Model, all_solvers: List[str], loopless: str
) -> None:
    """Test find_blocked_reactions() with loopless constraints."""
    model.solver = all_solvers
    regular_result = find_blocked_reactions(model, processes=1)
    loopless_result = find_blocked_reactions(model, loopless=loopless, processes=1)

    assert set(loopless_result) == set(regular_result)
    assert set(loopless_result.forward_blocked) == set(regular_result.forward_blocked)
    assert set(loopless_result.reverse_blocked) == set(regular_result.reverse_blocked)


@pytest.mark.parametrize("loopless", ["fastSNP", "potentials"])
def test_find_blocked_reactions_with_loopless_ll_test_model(
    ll_test_model: Model, loopless: str
) -> None:
    """Test find_blocked_reactions() with loopless constraints on ll_test_model."""
    loopless_result = find_blocked_reactions(
        ll_test_model, loopless=loopless, processes=1
    )

    assert loopless_result == ["v3"]
    assert loopless_result.forward_blocked == ["v3"]
    assert set(loopless_result.reverse_blocked) == {"EX_A", "DM_C", "v1", "v2", "v3"}
