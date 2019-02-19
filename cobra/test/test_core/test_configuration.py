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


def test_default_tolerances():
    """Verify the default solver tolerances."""
    config = Configuration()
    config.solver = "glpk"
    assert config.tolerances["feasibility"] == 1e-07
    assert config.tolerances["optimality"] == 1e-07
    assert config.tolerances["integrality"] == 1e-05


def test_tolerances():
    """Test assignment of solver tolerances."""
    config = Configuration()
    config.tolerances = {
        "feasibility": 1e-06,
        "optimality": 1e-06,
        "integrality": 1e-06,
    }
    assert config.tolerances["feasibility"] == 1e-06
    assert config.tolerances["optimality"] == 1e-06
    assert config.tolerances["integrality"] == 1e-06
    # Restore default tolerances
    config.tolerances = {
        "feasibility": 1e-07,
        "optimality": 1e-07,
        "integrality": 1e-05,
    }
