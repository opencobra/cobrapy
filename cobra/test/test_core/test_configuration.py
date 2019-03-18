# -*- coding: utf-8 -*-

"""Test functions of configuration.py"""

from __future__ import absolute_import

from cobra.core import Configuration, Model
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


def test_default_tolerance(model):
    """Verify the default solver tolerance."""
    config = Configuration()
    config.solver = "glpk"
    assert config.tolerance == 1e-07
    # Test the consistency between cobra.core.Configuration.tolerance and
    # cobra.core.Model.tolerance
    assert config.tolerance == model.tolerance


def test_toy_model_tolerance_with_different_default():
    """Verify that different default tolerance is respected by Model."""
    config = Configuration()
    config.tolerance = 1e-05

    toy_model = Model(name="toy model")
    assert toy_model.tolerance == 1e-05


def test_tolerance_assignment(model):
    """Test assignment of solver tolerance."""
    model.tolerance = 1e-06
    assert model.tolerance == 1e-06
