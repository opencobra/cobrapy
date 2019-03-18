# -*- coding: utf-8 -*-

"""Define the global configuration."""

from __future__ import absolute_import

import logging
import types
from multiprocessing import cpu_count
from warnings import warn

from six import string_types, with_metaclass

from cobra.core.singleton import Singleton
from cobra.exceptions import SolverNotFound
from cobra.util.solver import interface_to_str, solvers


__all__ = ("Configuration",)


LOGGER = logging.getLogger(__name__)


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
    tolerance: float
        The default tolerance for the solver being used (default 1E-09).
    lower_bound : float
        The standard lower bound for reversible reactions (default -1000).
    upper_bound : float
        The standard upper bound for all reactions (default 1000).
    bounds : tuple of floats
        The default reaction bounds for newly created reactions. The bounds
        are in the form of lower_bound, upper_bound (default -1000.0, 1000.0).
    processes : int
        A default number of processes to use where multiprocessing is
        possible. The default number corresponds to the number of available
        cores (hyperthreads).

    """

    def __init__(self):
        self._solver = None
        self.tolerance = 1E-07
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

        try:
            self.processes = cpu_count()
        except NotImplementedError:
            LOGGER.warning(
                "The number of cores could not be detected - assuming 1.")
            self.processes = 1

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
        return """
        solver: {solver}
        solver tolerance: {tolerance}
        lower_bound: {lower_bound}
        upper_bound: {upper_bound}
        processes: {processes}""".format(
            solver=interface_to_str(self.solver),
            tolerance=self.tolerance,
            lower_bound=self.lower_bound,
            upper_bound=self.upper_bound,
            processes=self.processes,
        )

    def _repr_html_(self):
        return """
        <table>
            <tr>
                <td><strong>Solver</strong></td>
                <td>{solver}</td>
            </tr>
            <tr>
                <td><strong>Solver tolerance</strong></td>
                <td>{tolerance}</td>
            </tr>
            <tr>
                <td><strong>Lower bound</strong></td>
                <td>{lower_bound}</td>
            </tr>
            <tr>
                <td><strong>Upper bound</strong></td>
                <td>{upper_bound}</td>
            </tr>
            <tr>
                <td><strong>Processes</strong></td>
                <td>{processes}</td>
            </tr>
        </table>""".format(
            solver=interface_to_str(self.solver),
            tolerance=self.tolerance,
            lower_bound=self.lower_bound,
            upper_bound=self.upper_bound,
            processes=self.processes,
        )


class Configuration(with_metaclass(Singleton, BaseConfiguration)):
    """Define the configuration to be singleton based."""

    pass
