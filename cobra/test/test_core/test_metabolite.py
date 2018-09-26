# -*- coding: utf-8 -*-

"""Test functions of metabolite.py"""

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


def test_repr_html_(model):
    assert '<table>' in model.metabolites.h2o_c._repr_html_()
