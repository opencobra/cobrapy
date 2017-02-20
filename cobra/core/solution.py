# -*- coding: utf-8 -*-

"""Provide unified interfaces to optimization solutions."""

from __future__ import absolute_import

import logging
from builtins import dict, object, super
from collections import OrderedDict
from warnings import warn

from cobra.exceptions import UndefinedSolution

__all__ = ("Solution",)

LOGGER = logging.getLogger(__name__)


class Solution(object):
    """
    The base for a unified interface to a `cobra.Model` optimization solution.

    Notes
    -----
    Solution values are lazily evaluated on a per request basis and thus
    accessing all `fluxes` and the like is rather inefficient. The lazy
    evaluation also means that all accessed values are stored for later use but
    those values that were not accessed may be invalidated by further
    optimizations of the same model.

    Attributes
    ----------
    model : cobra.Model
        The model used for finding a solution.
    objective_value : float
        The (optimal) value for the objective function.
    status : str
        The solver status related to the solution.
    fluxes : OrderedDict
        An ordered dictionary of fluxes.
    reduced_costs : OrderedDict
        An ordered dictionary of reduced costs.
    shadow_prices : OrderedDict
        An ordered dictionary of shadow prices.

    Deprecated Attributes
    ---------------------
    f : float
        Use `objective_value` instead.
    x : list
        Use `fluxes.values()` instead.
    x_dict : OrderedDict
        Use `fluxes` instead.
    y : list
        Use `reduced_costs.values()` instead.
    y_dict : OrderedDict
        Use `reduced_costs` instead.
    """

    def __new__(cls, *args, **kwargs):
        """Create either LegacySolution or Solution type."""
        # Prevent a circular import between Solution and Model.
        from cobra.core.model import Model
        # this is a cobrapy compatibility hack
        if len(args) == 1 and not isinstance(args[0], Model):
            cobrapy_solution = super(Solution, cls).__new__(LegacySolution)
            cobrapy_solution.__init__(*args, **kwargs)
            return cobrapy_solution
        else:
            return super(Solution, cls).__new__(cls)

    def __init__(self, model, **kwargs):
        """
        Initialize a unified solution interface from a model.

        Parameters
        ----------
        model : cobra.Model
            The model used for finding the solution.
        """
        super(Solution, self).__init__(**kwargs)
        self._model = model
        self._objective_value = None
        self._status = None
        self._fluxes = None
        self._reduced_costs = None
        self._shadow_prices = None
        self._store = dict()

    def __repr__(self):
        """String representation of the solution instance."""
        if self.objective_value is None:
            return "<Solution {0:r} at 0x{1:x}>".format(self.status, id(self))
        return "<Solution {0:.3g} at 0x{1:x}>".format(self.objective_value,
                                                      id(self))

    def __dir__(self):
        """Hide deprecated attributes and methods from the public interface."""
        fields = sorted(dir(type(self)) + list(self.__dict__))
        fields.remove('f')
        fields.remove('x')
        fields.remove('y')
        fields.remove('x_dict')
        fields.remove('y_dict')
        return fields

    def __getitem__(self, reaction_id):
        """
        Return the flux of a reaction.

        Parameters
        ----------
        reaction : str
            A model reaction ID.
        """
        flux = self._store.get(reaction_id)
        if flux is not None:
            return flux
        self._is_current()
        reaction = self._model.reactions.get_by_id(reaction_id)
        self._store[reaction_id] = flux = reaction.forward_variable.primal -\
            reaction.reverse_variable.primal
        return flux

    get_primal_by_id = __getitem__

    def _is_current(self):
        """
        Ensure that the solution is current.

        Raises
        ------
        UndefinedSolution
            If the solution has become invalid due to re-optimization of the
            underlying model.
        """
        if self is not self._model.solution:
            raise UndefinedSolution(
                "The solution {0:s} was invalidated by a more recent"
                " optimization solution {1:s}. Previously accesssed values are"
                " still available.".format(repr(self),
                                           repr(self._model.solution))
            )

    @property
    def objective_value(self):
        """Access the objective value."""
        if self._objective_value is not None:
            return self._objective_value
        self._is_current()
        self._objective_value = self._model.solver.objective.value
        return self._objective_value

    @property
    def status(self):
        """Access the solver status after optimization."""
        if self._status is not None:
            return self._status
        self._is_current()
        self._status = self._model.solver.status
        return self._status

    @property
    def fluxes(self):
        """
        Get the map of reaction IDs to fluxes.

        Warning
        -------
        Accessing all fluxes in this way is not recommended since it defeats
        the purpose of lazy evaluation.

        Returns
        -------
        OrderedDict
            All fluxes in the model as an ordered dictionary keyed by
            reaction ID.
        """
        if self._fluxes is not None:
            return self._fluxes
        self._is_current()
        self._fluxes = OrderedDict(
            (rxn.id, rxn.forward_variable.primal - rxn.reverse_variable.primal)
            for rxn in self._model.reactions)
        return self._fluxes

    @property
    def reduced_costs(self):
        """
        Get the map of reaction IDs to reduced costs.

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
        if self._reduced_costs is not None:
            return self._reduced_costs
        self._is_current()
        self._reduced_costs = OrderedDict(
            (rxn.id, rxn.forward_variable.dual - rxn.reverse_variable.dual)
            for rxn in self._model.reactions)
        return self._reduced_costs

    @property
    def shadow_prices(self):
        """
        Get the map of reaction IDs to shadow prices.

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
        if self._shadow_prices is not None:
            return self._shadow_prices
        self._is_current()
        self._shadow_prices = self._model.solver.shadow_prices
        return self._shadow_prices

    @property
    def f(self):
        """Deprecated property for getting the objective value."""
        warn("use solution.objective_value instead", DeprecationWarning)
        return self.objective_value

    @property
    def x_dict(self):
        """Deprecated property for getting fluxes."""
        warn("use solution.fluxes instead", DeprecationWarning)
        return self.fluxes

    @x_dict.setter
    def x_dict(self, fluxes):
        """Deprecated property for setting fluxes."""
        warn("let Model create a solution instance, don't update yourself",
             DeprecationWarning)
        self._fluxes = fluxes

    @property
    def x(self):
        """Deprecated property for getting flux values."""
        warn("use solution.fluxes.values() instead", DeprecationWarning)
        return self.fluxes.values()

    @property
    def y_dict(self):
        """Deprecated property for getting reduced costs."""
        warn("use solution.reduced_costs instead", DeprecationWarning)
        return self.reduced_costs

    @y_dict.setter
    def y_dict(self, costs):
        """Deprecated property for setting reduced costs."""
        warn("let Model create a solution instance, don't update yourself",
             DeprecationWarning)
        self._reduced_costs = costs

    @property
    def y(self):
        """Deprecated property for getting reduced cost values."""
        warn("use solution.reduced_costs.values() instead", DeprecationWarning)
        return self.reduced_costs.values()


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

    def __init__(self, f, x=None, x_dict=None, y=None, y_dict=None,
                 solver=None, the_time=0, status='NA', **kwargs):
        """
        Initialize a legacy interface to a solution from an objective value.

        Parameters
        ----------
        f : float
            Objective value.
        solver : str, optional
            A string indicating which solver package was used.
        x : iterable, optional
            List or Array of the fluxes (primal values).
        x_dict : dict, optional
            A dictionary of reaction IDs that maps to the respective primal
            values.
        y : iterable, optional
            List or Array of the dual values.
        y_dict : dict, optional
            A dictionary of reaction IDs that maps to the respective dual
            values.
        the_time : int, optional
        status : str, optional

        .. warning :: deprecated
        """
        super(LegacySolution, self).__init__(**kwargs)
        self.solver = solver
        self.f = f
        self.x = x
        self.x_dict = x_dict
        self.status = status
        self.y = y
        self.y_dict = y_dict

    def __repr__(self):
        """String representation of the solution instance."""
        if self.objective_value is None:
            return "<LegacySolution {0:r} at 0x{1:x}>".format(
                self.status, id(self))
        return "<LegacySolution {0:.3g} at 0x{1:x}>".format(
            self.objective_value, id(self))

    def __getitem__(self, reaction_id):
        """
        Return the flux of a reaction.

        Parameters
        ----------
        reaction_id : str
            A reaction ID.
        """
        return self.x_dict[reaction_id]

    def dress_results(self, model):
        """
        Method could be intended as a decorator.

        .. warning :: deprecated
        """
        warn("unnecessary to call this deprecated function",
             DeprecationWarning)
