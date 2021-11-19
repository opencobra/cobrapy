# -*- coding: utf-8 -*-

"""Test functions of gene.py"""

from cobra.core.gene import GPR


def test_repr_html_(model):
    assert "<table>" in model.genes[0]._repr_html_()


def test_gpr():
    gpr1 = GPR()
    gpr1.from_string(string_gpr="")
    assert isinstance(gpr1.body, list)
    assert len(gpr1.body) == 0


def test_gpr_functions():
    gpr1 = GPR()
    # Completely empty GPRs will return True
    assert gpr1.eval()
    assert gpr1.to_string() == ""
    gpr2 = gpr1.copy()
