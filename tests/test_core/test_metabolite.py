"""Test functions of metabolite.py."""

from typing import TYPE_CHECKING, Tuple

import pytest

from cobra.core import Metabolite


if TYPE_CHECKING:
    from cobra import Model, Solution


def test_metabolite_formula() -> None:
    """Test metabolite formula correctly leads to element dictionary."""
    met = Metabolite("water")
    met.formula = "H2O"
    assert met.elements == {"H": 2, "O": 1}
    assert met.formula_weight == 18.01528


def test_formula_element_setting(model: "Model") -> None:
    """Test that formula and elements set each other when one is set."""
    met = model.metabolites[1]
    orig_formula = str(met.formula)
    orig_elements = dict(met.elements)
    met.formula = ""
    assert met.elements == {}
    met.elements = orig_elements
    assert met.formula == orig_formula


def test_set_id(solved_model: Tuple["Solution", "Model"]) -> None:
    """Test that setting id leads to change in the model metabolite dictlist."""
    solution, model = solved_model
    met = Metabolite("test")
    with pytest.raises(TypeError):
        met.id = 1
    model.add_metabolites([met])
    with pytest.raises(ValueError):
        met.id = "g6p_c"
    met.id = "test2"
    assert "test2" in model.metabolites
    assert "test" not in model.metabolites


def test_remove_from_model(solved_model: Tuple["Solution", "Model"]) -> None:
    """Test that removing a model removes it from relevant model dictlists."""
    solution, model = solved_model
    met = model.metabolites.get_by_id("g6p_c")
    met.remove_from_model()
    assert not (met.id in model.metabolites)
    assert not (met.id in model.constraints)


def test_repr_html_(model: "Model") -> None:
    """Test HTML represenation is correct for a metabolite."""
    assert "<table>" in model.metabolites.h2o_c._repr_html_()
