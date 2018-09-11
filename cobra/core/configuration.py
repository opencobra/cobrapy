# -*- coding: utf-8 -*-

"""Contains class for configuration."""

from __future__ import absolute_import

from cobra.core import Model, Reaction


class Configuration(object):
    """
    A class for setting and storing global configuration values.

    This object sets default values for the modifiable attributes like
    default solver, reaction bounds etc.

    Attributes
    ----------
    solver: str, optional (default "glpk")
        The default solver. The choices are the ones provided
        by optlang: "glpk", "cplex" and "gurobi".
    bounds: tuple of floats, optional (default (0.0, 1000.0))
        The default reaction bounds. The reaction bounds are in the form
        of (lower_bound, upper_bound).

    """
    __instance = None

    def __new__(cls, *args, **kwargs):
        if cls.__instance is None:
            cls.__instance = object.__new__(cls)
        return cls.__instance


    def __init__(self, solver="glpk", bounds=(0.0, 1000.0)):
        self._solver = solver
        self._bounds = bounds


    def set(self):
        """Bind the values to the respective objects."""
        # TODO: complete the binding
        # Model.solver = self._solver
        # Reaction.bounds = self._bounds


    @property
    def solver(self):
        return self._solver


    @property
    def default_bounds(self):
        return self._bounds


    def __repr__(self):
        return "<Configuration at 0x{}>".format(id(self))
