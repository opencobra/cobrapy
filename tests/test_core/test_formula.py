"""Test functions of formula.py ."""

import pytest

from cobra.core.formula import Formula


def test_formula_init() -> None:
    """Test initialization."""
    f = Formula("H2O")
    assert "Formula H2O" in repr(f)
    f = Formula("H2O", name="Water")
    assert f.name == "Water"


@pytest.mark.parametrize(
    ["formula", "weight"], [["H2O", 18.01528], ["C6H12O6", 180.15588], ["NO3", 62.0049]]
)
def test_formula_weight(formula, weight) -> None:
    """Test molecular weight calculation."""
    assert Formula(formula).weight == pytest.approx(weight)


def test_formula_wrong() -> None:
    """Test incorrect formula elements."""
    with pytest.warns(UserWarning):
        w = Formula("NOX").weight
        assert w is None
