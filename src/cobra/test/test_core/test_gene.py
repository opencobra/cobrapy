# -*- coding: utf-8 -*-

"""Test functions of gene.py"""

from cobra.core.gene import gpr_eq


def test_repr_html_(model):
    assert "<table>" in model.genes[0]._repr_html_()


def test_gpr_eq():
    assert gpr_eq("a and b", "b and a")
    assert gpr_eq("a AND b", "a AND b")
    assert gpr_eq("", "")
    assert gpr_eq("(a and b) or (a and c)", "a and (b or c)")
    assert gpr_eq("(a and b) or (a and c)", "(a and c) or (a and b)")
