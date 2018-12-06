# -*- coding: utf-8 -*-

"""Test functions of metabolite.py"""

import pytest

from cobra.core import Metabolite


def test_metabolite_formula():
    met = Metabolite("water")
    met.formula = "H2O"
    assert met.elements == {"H": 2, "O": 1}
    assert met.formula_weight == 18.01528


def test_formula_element_setting(model):
    met = model.metabolites[1]
    orig_formula = str(met.formula)
    orig_elements = dict(met.elements)
    met.formula = ''
    assert met.elements == {}
    met.elements = orig_elements
    assert met.formula == orig_formula


def test_set_id(solved_model):
    solution, model = solved_model
    met = Metabolite("test")
    with pytest.raises(TypeError):
        setattr(met, 'id', 1)
    model.add_metabolites([met])
    with pytest.raises(ValueError):
        setattr(met, "id", 'g6p_c')
    met.id = "test2"
    assert "test2" in model.metabolites
    assert "test" not in model.metabolites


def test_remove_from_model(solved_model):
    solution, model = solved_model
    met = model.metabolites.get_by_id("g6p_c")
    met.remove_from_model()
    assert not (met.id in model.metabolites)
    assert not (met.id in model.constraints)


def test_repr_html_(model):
    assert '<table>' in model.metabolites.h2o_c._repr_html_()
