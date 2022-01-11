# -*- coding: utf-8 -*-

"""Test functions of gene.py"""

from cobra.core.gene import GPR


def test_repr_html_(model):
    assert "<table>" in model.genes[0]._repr_html_()
