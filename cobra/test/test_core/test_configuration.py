# -*- coding: utf-8 -*-

"""Test functions of configuration.py"""

from __future__ import absolute_import

from cobra.core import Configuration
from cobra.util.solver import interface_to_str


def test_default_bounds():
    """Verify the default bounds."""
    config = Configuration()
    assert config.bounds == (-1000.0, 1000.0)


def test_bounds():
    """Test changing bounds."""
    config = Configuration()
    config.bounds = 100.0, 10000.0
    assert config.bounds == (100.0, 10000.0)
    # Restore default values.
    config.bounds = -1000.0, 1000.0


def test_solver():
    """Test assignment of different solvers."""
    config = Configuration()
    config.solver = "glpk"
    assert interface_to_str(config.solver) == "glpk"
    config.solver = "glpk_exact"
    assert interface_to_str(config.solver) == "glpk_exact"
    # Restore default solver.
    config.solver = "glpk"
