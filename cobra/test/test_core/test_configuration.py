# -*- coding: utf-8 -*-

"""Test functions of configuration.py"""

from __future__ import absolute_import

from cobra.core import Configuration
from cobra.util.solver import interface_to_str


def test_bounds():
    config = Configuration()
    # Check default bounds
    assert config.bounds == (0.0, 1000.0)
    # Change bounds
    config.bounds = (100.0, 10000.0)
    # Check new bounds
    assert config.bounds == (100.0, 10000.0)


def test_solver():
    config = Configuration()
    # Check default solver
    assert interface_to_str(config.solver) == "glpk"
    # Change solver
    config.solver = "glpk_exact"
    # Check new solver
    assert interface_to_str(config.solver) == "glpk_exact"
    # Change it back for further tests
    config.bounds = (0.0, 100.0)
    config.solver = "glpk"
