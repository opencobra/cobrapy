"""Provide a global configuration object."""


import logging
import pathlib
import types
from numbers import Number
from os import cpu_count
from textwrap import dedent
from typing import Optional, Tuple, Union

import appdirs

from ..exceptions import SolverNotFound
from ..util.solver import interface_to_str
from ..util.solver import solvers as SOLVERS
from .singleton import Singleton


__all__ = ("Configuration",)


logger = logging.getLogger(__name__)


class Configuration(metaclass=Singleton):
    """
    Define a global configuration object.

    The attributes of this singleton object are used as default values by cobra
    functions.

    Attributes
    ----------
    solver : {"glpk", "cplex", "gurobi", "glpk_exact"}
        The default solver for new models. The solver choices are the ones
        provided by `optlang` and depend on solvers installed in your environment.
    tolerance : float, optional
        The default tolerance for the solver being used (default 1E-07).
    lower_bound : float, optional
        The standard lower bound for reversible reactions (default -1000).
    upper_bound : float, optional
        The standard upper bound for all reactions (default 1000).
    bounds : tuple of floats
        The default reaction bounds for newly created reactions. The bounds
        are in the form of lower_bound, upper_bound (default -1000.0, 1000.0).
    processes : int > 0
        A default number of processes to use where multiprocessing is
        possible. The default number corresponds to the number of available
        cores (hyperthreads) minus one.
    cache_directory : pathlib.Path or str, optional
        A path where the model cache should reside if caching is desired. The
        default directory depends on the operating system.
    max_cache_size : int, optional
        The allowed maximum size of the model cache in bytes (default 1 GB).
    cache_expiration : int, optional
        The expiration time in seconds for the model cache if any (default None).

    """

    def __init__(self, **kwargs) -> None:
        """Initialize the configuration with its default attribute values."""
        super().__init__(**kwargs)
        self._solver = None
        self.tolerance = 1e-07
        self.lower_bound = None
        self.upper_bound = None
        self.processes = None
        self._cache_directory = None
        # Set the cache size to a maximum of 100 MB.
        self.max_cache_size = 100 * (1024**2)
        self.cache_expiration = None

        self.bounds = -1000.0, 1000.0
        self._set_default_solver()
        self._set_default_processes()
        self._set_default_cache_directory()

    def _set_default_solver(self) -> None:
        """Set the default solver from a preferred order."""
        for name in ["gurobi", "cplex", "glpk"]:
            try:
                self.solver = name
            except SolverNotFound:
                continue
            else:
                break

    def _set_default_processes(self) -> None:
        """Set the default number of processes."""
        self.processes = cpu_count()
        if self.processes is None:
            logger.warning("The number of cores could not be detected - assuming one.")
            self.processes = 1
        if self.processes > 1:
            self.processes -= 1

    def _set_default_cache_directory(self) -> None:
        """Set the platform-dependent default cache directory."""
        self.cache_directory = pathlib.Path(
            appdirs.user_cache_dir(appname="cobrapy", appauthor="opencobra")
        )

    @property
    def solver(self) -> types.ModuleType:
        """Return the optlang solver interface."""
        return self._solver

    @solver.setter
    def solver(self, value) -> None:
        """Set the optlang solver interface."""
        not_valid_interface = SolverNotFound(
            f"'{value}' is not a valid solver interface. "
            f" Please pick one from {', '.join(SOLVERS)}."
        )
        if isinstance(value, str):
            if value not in SOLVERS:
                raise not_valid_interface
            interface = SOLVERS[value]
        elif isinstance(value, types.ModuleType) and hasattr(value, "Model"):
            interface = value
        else:
            raise not_valid_interface
        self._solver = interface

    @property
    def bounds(self) -> Tuple[Optional[Number], Optional[Number]]:
        """Return the lower, upper reaction bound pair.

        Returns
        -------
        tuple of number and number or None and None
            The lower and upper bounds for new reactions.

        """
        return self.lower_bound, self.upper_bound

    @bounds.setter
    def bounds(self, bounds: Tuple[Optional[Number], Optional[Number]]) -> None:
        """Set the lower, upper reaction bound pair.

        Parameters
        ----------
        bounds : tuple of number and number or None and None
            The lower and upper bounds for new reactions.

        """
        if None not in bounds:
            assert bounds[0] <= bounds[1]
        self.lower_bound = bounds[0]
        self.upper_bound = bounds[1]

    @property
    def cache_directory(self) -> pathlib.Path:
        """Return the model cache directory."""
        return self._cache_directory

    @cache_directory.setter
    def cache_directory(self, path: Union[pathlib.Path, str]) -> None:
        """
        Set the model cache directory.

        The directory path is created if it doesn't exist yet.

        Parameters
        ----------
        path : pathlib.Path or str
            The path to the cache directory.

        """
        self._cache_directory = pathlib.Path(path)
        if not self._cache_directory.is_dir():
            logger.debug(f"Creating cache directory '{str(self._cache_directory)}'.")
            self._cache_directory.mkdir(parents=True)

    def __repr__(self) -> str:
        """Return a string representation of the current configuration values."""
        return dedent(
            f"""
            solver: {interface_to_str(self.solver)}
            tolerance: {self.tolerance}
            lower_bound: {self.lower_bound}
            upper_bound: {self.upper_bound}
            processes: {self.processes}
            cache_directory: {self.cache_directory}
            max_cache_size: {self.max_cache_size}
            cache_expiration: {self.cache_expiration}
            """
        )

    def _repr_html_(self) -> str:
        """
        Return a rich HTML representation of the current configuration values.

        Notes
        -----
        This special method is used automatically in Jupyter notebooks to display a
        result from a cell.

        """
        return dedent(
            f"""
            <table>
              <thead>
                <tr>
                  <td><strong>Attribute</strong></td>
                  <td><strong>Description</strong></td>
                  <td><strong>Value</strong></td>
                </tr>
              </thead>
              <tbody>
                <tr>
                  <td><pre>solver</pre></td>
                  <td>Mathematical optimization solver</td>
                  <td>{interface_to_str(self.solver)}</td>
                </tr>
                <tr>
                    <td><pre>tolerance</pre></td>
                    <td>General solver tolerance (feasibility, integrality, etc.)</td>
                    <td>{self.tolerance}</td>
                </tr>
                <tr>
                    <td><pre>lower_bound</pre></td>
                    <td>Default reaction lower bound</td>
                    <td>{self.lower_bound}</td>
                </tr>
                <tr>
                    <td><pre>upper_bound</pre></td>
                    <td>Default reaction upper bound</td>
                    <td>{self.upper_bound}</td>
                </tr>
                <tr>
                    <td><pre>processes</pre></td>
                    <td>Number of parallel processes</td>
                    <td>{self.processes}</td>
                </tr>
                <tr>
                    <td><pre>cache_directory</pre></td>
                    <td>Path for the model cache</td>
                    <td>{self.cache_directory}</td>
                </tr>
                <tr>
                    <td><pre>max_cache_size</pre></td>
                    <td>Maximum cache size in bytes</td>
                    <td>{self.max_cache_size}</td>
                </tr>
                <tr>
                    <td><pre>cache_expiration</pre></td>
                    <td>Model cache expiration time in seconds (if any)</td>
                    <td>{self.cache_expiration}</td>
                </tr>
              </tbody>
            </table>
            """
        )
