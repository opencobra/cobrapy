"""Test functionalities of flux analysis helper functions."""

import pytest

from cobra.core import Model
from cobra.flux_analysis.helpers import normalize_cutoff


def test_normalize_cutoff(model: Model) -> None:
    """Test normalize cutoff."""
    cutoff = normalize_cutoff(model)
    assert cutoff == 1e-7


def test_normalize_cutoff_with_specified_cutoff_above_default(
    model: Model,
) -> None:
    """Test normalize cutoff with specified cutoff greater than default."""
    cutoff = normalize_cutoff(model, 1e-3)
    assert cutoff == 1e-3


def test_normalize_cutoff_with_specified_cutoff_below_default(
    model: Model,
) -> None:
    """Test normalize cutoff with specified cutoff less than default."""
    with pytest.raises(ValueError):
        normalize_cutoff(model, 1e-10)
