# -*- coding: utf-8 -*-

"""Provide unified interfaces to optimization solutions."""

from __future__ import absolute_import

import logging
from builtins import object, super
from warnings import warn

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
    objective_value : float
        The (optimal) value for the objective function.
    status : str
        The solver status related to the solution.
    fluxes : pandas.Series
        Contains the reaction fluxes (primal values of variables).
    reduced_costs : pandas.Series
        Contains reduced costs (dual values of variables).

    Deprecated Attributes
    ---------------------
    f : float
        Use `objective_value` instead.
    x : list
        Use `fluxes.values` instead.
    x_dict : pandas.Series
        Use `fluxes` instead.
    y : list
        Use `reduced_costs.values` instead.
    y_dict : pandas.Series
        Use `reduced_costs` instead.
    """

    def __init__(self, reactions, objective_value, status, fluxes,
                 reduced_costs=None, shadow_prices=None, **kwargs):
        """
        Initialize a unified solution interface from a model.

        Parameters
        ----------
        reactions : iterable
            A list of `cobra.Reaction` objects for which the solution is
            retrieved.
        objective_value : float
            The (optimal) value for the objective function.
        status : str
            The solver status related to the solution.
        fluxes : dict-like
            A dict-like container that is indexable by `cobra.Reaction` objects or
            their IDs, for example, `OrderedDict`, `cobra.DictList`,
            `optlang.Container`, or eventually `pandas.Series`. Contains reaction
            fluxes (primal values).
        reduced_costs : dict-like
            As with fluxes but contains reduced costs (dual values).
        shadow_prices : dict-like
            As with fluxes but contains shadow prices.
        """
        super(Solution, self).__init__(**kwargs)
        self.reactions = reactions
        self.objective_value = objective_value
        self.status = status
        self.fluxes = fluxes
        self.reduced_costs = reduced_costs
        self.shadow_prices = shadow_prices

    def __repr__(self):
        """String representation of the solution instance."""
        if self.status != "optimal":
            return "<Solution {0:r} at 0x{1:x}>".format(self.status, id(self))
        return "<Solution {0:g} at 0x{1:x}>".format(self.objective_value,
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
        return self.fluxes[reaction_id]

    get_primal_by_id = __getitem__

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
        warn("let Model.optimize create a solution instance, don't update yourself",
             DeprecationWarning)
        self.fluxes = fluxes

    @property
    def x(self):
        """Deprecated property for getting flux values."""
        warn("use solution.fluxes.values() instead", DeprecationWarning)
        return self.fluxes.values

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
        self.reduced_costs = costs

    @property
    def y(self):
        """Deprecated property for getting reduced cost values."""
        warn("use solution.reduced_costs.values() instead", DeprecationWarning)
        return self.reduced_costs.values


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

    Warning
    -------
    The LegacySolution class and its interface is deprecated.
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
        if self.status != "optimal":
            return "<LegacySolution {0:r} at 0x{1:x}>".format(
                self.status, id(self))
        return "<LegacySolution {0:g} at 0x{1:x}>".format(
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
