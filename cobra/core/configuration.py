# -*- coding: utf-8 -*-

"""Define the global configuration."""

from __future__ import absolute_import

import types

from six import with_metaclass, string_types

from cobra.exceptions import SolverNotFound
from cobra.core.singleton import Singleton
from cobra.util.solver import interface_to_str, solvers


__all__ = ("Configuration",)


class BaseConfiguration(object):
    """
    Define global configuration value that will be honored by cobra functions.

    This object sets default values for the modifiable attributes like
    default solver, reaction bounds etc.

    Attributes
    ----------
    solver : {"glpk", "cplex", "gurobi"}
        The default solver for new models. The solver choices are the ones
        provided by `optlang` and solvers installed in your environment.
    lower_bound : float
        The standard lower bound for reversible reactions (default -1000).
    upper_bound : float
        The standard upper bound for all reactions (default 1000).
    bounds : tuple of floats
        The default reaction bounds for newly created reactions. The bounds
        are in the form of lower_bound, upper_bound (default -1000.0, 1000.0).

    """

    def __init__(self):
        self._solver = None
        self.lower_bound = None
        self.upper_bound = None
        # Set the default solver from a preferred order.
        for name in ["gurobi", "cplex", "glpk"]:
            try:
                self.solver = name
            except SolverNotFound:
                continue
            else:
                break
        self.bounds = -1000.0, 1000.0

    @property
    def solver(self):
        return self._solver

    @solver.setter
    def solver(self, value):
        not_valid_interface = SolverNotFound(
            "'{}' is not a valid solver interface. Pick one from {}.".format(
                value, ", ".join(list(solvers))))
        if isinstance(value, string_types):
            if value not in solvers:
                raise not_valid_interface
            interface = solvers[value]
        elif isinstance(value, types.ModuleType) and hasattr(value, 'Model'):
            interface = value
        else:
            raise not_valid_interface
        self._solver = interface

    @property
    def bounds(self):
        return self.lower_bound, self.upper_bound

    @bounds.setter
    def bounds(self, bounds):
        # TODO: We should consider allowing `None` for free bounds.
        assert bounds[0] <= bounds[1]
        self.lower_bound = bounds[0]
        self.upper_bound = bounds[1]

    def __repr__(self):
        return "solver: {}\nlower_bound: {}\nupper_bound: {}".format(
            interface_to_str(self.solver), self.lower_bound, self.upper_bound)


class Configuration(with_metaclass(Singleton, BaseConfiguration)):
    """Define the configuration to be singleton based."""

    pass
