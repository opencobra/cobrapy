# -*- coding: utf-8 -*-

"""
Provide classes for unified solver solution interfaces.

`cobra.solution`
`cobra.LazySolution`
`cobra.LegacySolution`
"""

from __future__ import absolute_import

import time
import datetime
import logging

from collections import OrderedDict
from warnings import warn

import cobra

try:
    import pandas
except ImportError:
    pandas = None

from cobra.exceptions import UndefinedSolution

logger = logging.getLogger(__name__)


class SolutionBase(object):
    """
    The base for a unified interface to a `cobra.Model` optimization solution.

    Attributes
    ----------
    model : cobra.Model
        The model used for finding a solution.
    objective_value : float
        The (optimal) value for the objective function.
    """

    def __new__(cls, *args, **kwargs):
        """Create either legacy or SolutionBase."""
        # this is a cobrapy compatibility hack
        if len(args) == 1 and not isinstance(args[0], cobra.core.Model):
            cobrapy_solution = super(SolutionBase, cls).__new__(LegacySolution)
            cobrapy_solution.__init__(*args, **kwargs)
            return cobrapy_solution
        else:
            return super(SolutionBase, cls).__new__(cls)

    def __init__(self, model, *args, **kwargs):
        """
        Initialize a base solution from a model.

        Parameters
        ----------
        model : cobra.Model
            The model used for finding a solution.
        """
        self.f = None
        self.model = model
        self._x = None
        self._y = None
        self._x_dict = None
        self._y_dict = None

    def get_primal_by_id(self, reaction_id):
        """Return a flux/primal value for a reaction.

        Parameters
        ----------
        reaction_id : str
            A reaction ID.
        """
        return self.x_dict[reaction_id]

    @property
    def x_dict(self):
        """Deprecated property for accessing fluxes."""
        warn("use solution.fluxes instead", DeprecationWarning)
        if self._x_dict is None:
            return self.fluxes
        else:
            return self._x_dict

    @x_dict.setter
    def x_dict(self, value):
        """Deprecated property for setting fluxes."""
        warn("not used", DeprecationWarning)
        self._x_dict = value

    @property
    def x(self):
        """Deprecated property for accessing flux values."""
        warn("use solution.fluxes.values instead", DeprecationWarning)
        if self._x is None:
            return self.fluxes.values()
        else:
            return self._x

    @x.setter
    def x(self, value):
        """Deprecated property for setting flux values."""
        warn("not used", DeprecationWarning)
        self._x = value

    @property
    def y_dict(self):
        """Deprecated property for accessing reduced costs."""
        warn("use solution.reduced_costs instead", DeprecationWarning)
        if self._y_dict is None:
            return self.reduced_costs
        else:
            return self._y_dict

    @y_dict.setter
    def y_dict(self, value):
        """Deprecated property for setting reduced costs."""
        self._y_dict = value

    @property
    def y(self):
        """Deprecated property for accessing reduced cost values."""
        warn("use solution.reduced_costs.values instead", DeprecationWarning)
        if self._y is None:
            return self.reduced_costs.values()
        else:
            return self._y

    @y.setter
    def y(self, value):
        """Deprecated property for setting reduced cost values."""
        warn("not used", DeprecationWarning)
        self._y = value

    @property
    def objective_value(self):
        """Return the objective value if it exists."""
        return self.f

    def __repr__(self):
        """String representation of the solution instance."""
        if self.f is None:
            return "<Solution {0:r} at 0x{1:x}>".format(self.status, id(self))
        return "<Solution {0:.3g} at 0x{1:x}>" % (self.f, id(self))


class Solution(SolutionBase):
    """
    A unified interface to a `cobra.Model` optimization solution.

    This is used to provide a single interface to results from different
    solvers that store their values in different ways.

    Attributes
    ----------
    fluxes : OrderedDict
        A dictionary of flux values.
    reduced_costs : OrderedDict
        A dictionary of reduced costs.
    shadow_prices : OrderedDict
    status : str


    Notes
    -----
    See also documentation for cobra.core.solution.solutionBase for a list of
    inherited attributes.
    """

    def __init__(self, model, *args, **kwargs):
        """
        Initialize a convenient solution interface from a model.

        Parameters
        ----------
        model : cobra.Model
            The model used for finding a solution.
        """
        super(Solution, self).__init__(model, *args, **kwargs)
        self.f = model.solver.objective.value
        self.fluxes = OrderedDict()
        self.shadow_prices = model.solver.shadow_prices
        self.reduced_costs = OrderedDict()
        self._primal_values = model.solver.primal_values
        self._reduced_values = model.solver.reduced_costs
        for reaction in model.reactions:
            self.fluxes[reaction.id] = (
                self._primal_values[reaction._get_forward_id()] -
                self._primal_values[reaction._get_reverse_id()]
            )
            self.reduced_costs[reaction.id] = (
                self._reduced_values[reaction._get_forward_id()] -
                self._reduced_values[reaction._get_reverse_id()]
            )

        self.status = model.solver.status
        self._reaction_ids = [r.id for r in self.model.reactions]
        self._metabolite_ids = [m.id for m in self.model.metabolites]

    def __dir__(self):
        """Hide deprecated attributes and methods from the public interface."""
        fields = sorted(dir(type(self)) + list(self.__dict__.keys()))
        fields.remove('x')
        fields.remove('y')
        fields.remove('x_dict')
        fields.remove('y_dict')
        return fields


class LazySolution(SolutionBase):
    """
    A lazily evaluated interface to a `cobra.Model` optimization solution.

    Instead of directly fetching results from the solver, this class only
    gets results when they are requested after checking that the model has
    not changed since the last optimization.

    Attributes
    ----------
    model : cobra.Model
        The model used for finding a solution.
    fluxes : OrderedDict
        A dictionary with flux values, populated upon first request.
    reduced_costs : OrderedDict
        A dictionary with the reduced costs for each reaction, populated
        upon first request.
    shadow_prices: OrderedDict
        A dictionary with the shadow_prices for each reaction, populated
        upon first request.

    Notes
    -----
    See also documentation for `cobra.core.solution.solutionBase` for an
    extensive list of inherited attributes.

    """

    def __init__(self, model, *args, **kwargs):
        """
        Initialize a lazily evaluated solution interface from a model.

        Parameters
        ----------
        model : cobra.Model
            The model used for finding a solution.
        """
        super(LazySolution, self).__init__(model, *args, **kwargs)
        if self.model._timestamp_last_optimization is not None:
            self._time_stamp = self.model._timestamp_last_optimization
        else:
            self._time_stamp = time.time()
        self._f = None
        self._primal_values = None
        self._reduced_values = None

    @property
    def data_frame(self):
        """Return flux values and reduced costs as a `pandas.DataFrame`."""
        if pandas is not None:
            return pandas.DataFrame(
                {'fluxes': pandas.Series(self.fluxes),
                 'reduced_costs': pandas.Series(self.reduced_costs)})
        else:
            warn("pandas not available")

    def _repr_html_(self):
        """Create an HTML representation of the solution, useful for Jupyter."""
        if pandas:
            return self.data_frame._repr_html_()
        else:
            warn("pandas not available")

    def _check_freshness(self):
        """
        Ensure that the solution is current.

        Raises
        ------
        UndefinedSolution
            If the solution has become invalid due to re-optimization of the
            underlying model.
        """
        # Assume that self.model._timestamp_last_optimization is not None since
        # otherwise there would be no solution.
        if self._time_stamp != self.model._timestamp_last_optimization:
            def timestamp_formatter(timestamp):
                datetime.datetime.fromtimestamp(timestamp).strftime(
                    "%Y-%m-%d %H:%M:%S:%f")

            raise UndefinedSolution(
                "The solution (captured around {0}) has become invalid as the "
                "model has been re-optimized recently ({1}).".format(
                    timestamp_formatter(self._time_stamp),
                    timestamp_formatter(self.model._timestamp_last_optimization)
                )
            )

    @property
    def status(self):
        """Access the solver status after optimization."""
        self._check_freshness()
        return self.model.solver.status

    @property
    def f(self):
        """Access the objective value."""
        self._check_freshness()
        if self._f is None:
            return self.model.solver.objective.value
        else:
            return self._f

    @f.setter
    def f(self, value):
        """Set the objective value."""
        self._f = value

    @property
    def fluxes(self):
        """
        Access the fluxes.

        Warning
        -------
        Accessing all flux values in this way is not recommended since it
        defeats the purpose of lazy evaluation.

        Returns
        -------
        OrderedDict
            All fluxes in the model as an ordered dictionary keyed by
            reaction ID.
        """
        self._check_freshness()
        primal_values = self.model.solver.primal_values

        fluxes = OrderedDict()
        for reaction in self.model.reactions:
            fluxes[reaction.id] = (
                primal_values[reaction._get_forward_id()] -
                primal_values[reaction._get_reverse_id()]
            )

        return fluxes

    @property
    def reduced_costs(self):
        """
        Access the reduced costs.

        Warning
        -------
        Accessing all reduced costs in this way is not recommended since it
        defeats the purpose of lazy evaluation.

        Returns
        -------
        OrderedDict
            All reduced costs in the model as an ordered dictionary keyed by
            reaction ID.
        """
        self._check_freshness()
        reduced_values = self.model.solver.reduced_costs

        reduced_costs = OrderedDict()
        for reaction in self.model.reactions:
            reduced_costs[reaction.id] = (
                reduced_values[reaction._get_forward_id()] -
                reduced_values[reaction._get_reverse_id()]
            )
        return reduced_costs

    @property
    def shadow_prices(self):
        """
        Access shadow prices.

        Warning
        -------
        Accessing all shadow prices in this way is not recommended since it
        defeats the purpose of lazy evaluation.

        Returns
        -------
        OrderedDict
            All shadow prices in the model as an ordered dictionary keyed by
            reaction ID.
        """
        self._check_freshness()
        return self.model.solver.shadow_prices

    def get_primal_by_id(self, reaction_id):
        """Return a flux/primal value for a reaction.

        Parameters
        ----------
        reaction_id : str
            A reaction ID.
        """
        self._check_freshness()
        return self.model.reactions.get_by_id(reaction_id).flux


class LegacySolution(object):
    """
    Legacy support for an interface to a `cobra.Model` optimization solution.

    Attributes
    ----------
    f : float
        The objective value
    solver : str
        A string indicating which solver package was used.
    x : iterable
        List or Array of the fluxes (primal values).
    x_dict : dict
        A dictionary of reaction IDs that maps to the respective primal values.
    y : iterable
        List or Array of the dual values.
    y_dict : dict
        A dictionary of reaction IDs that maps to the respective dual values.

    .. warning :: deprecated
    """

    def __init__(self, f, x=None,
                 x_dict=None, y=None, y_dict=None,
                 solver=None, the_time=0, status='NA'):
        """
        Initialize a legacy interface to a solution from an objective value.

        Parameters
        ----------
        f : float
            Objective value.

        .. warning :: deprecated
        """
        self.solver = solver
        self.f = f
        self.x = x
        self.x_dict = x_dict
        self.status = status
        self.y = y
        self.y_dict = y_dict

    def dress_results(self, model):
        """
        Method could be intended as a decorator.

        .. warning :: deprecated
        """
        warn("unnecessary to call this deprecated function",
             DeprecationWarning)
