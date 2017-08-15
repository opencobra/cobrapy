# -*- coding: utf-8 -*-

"""Provide unified interfaces to optimization solutions."""

from __future__ import absolute_import

import logging
from builtins import object, super
from warnings import warn

from numpy import empty, nan
from optlang.interface import OPTIMAL
from pandas import Series, DataFrame, option_context

from cobra.util.solver import check_solver_status

__all__ = ("Solution", "LegacySolution", "get_solution")

LOGGER = logging.getLogger(__name__)


class Solution(object):
    """
    A unified interface to a `cobra.Model` optimization solution.

    Notes
    -----
    Solution is meant to be constructed by `get_solution` please look at that
    function to fully understand the `Solution` class.

    Attributes
    ----------
    objective_value : float
        The (optimal) value for the objective function.
    status : str
        The solver status related to the solution.
    fluxes : pandas.Series
        Contains the reaction fluxes (primal values of variables).
    reduced_costs : pandas.Series
        Contains reaction reduced costs (dual values of variables).
    shadow_prices : pandas.Series
        Contains metabolite shadow prices (dual values of constraints).

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

    def __init__(self, objective_value, status, fluxes, reduced_costs=None,
                 shadow_prices=None, **kwargs):
        """
        Initialize a `Solution` from its components.

        Parameters
        ----------
        objective_value : float
            The (optimal) value for the objective function.
        status : str
            The solver status related to the solution.
        fluxes : pandas.Series
            Contains the reaction fluxes (primal values of variables).
        reduced_costs : pandas.Series
            Contains reaction reduced costs (dual values of variables).
        shadow_prices : pandas.Series
            Contains metabolite shadow prices (dual values of constraints).
        """
        super(Solution, self).__init__(**kwargs)
        self.objective_value = objective_value
        self.status = status
        self.fluxes = fluxes
        self.reduced_costs = reduced_costs
        self.shadow_prices = shadow_prices

    def __repr__(self):
        """String representation of the solution instance."""
        if self.status != OPTIMAL:
            return "<Solution {0:s} at 0x{1:x}>".format(self.status, id(self))
        return "<Solution {0:.3f} at 0x{1:x}>".format(self.objective_value,
                                                      id(self))

    def _repr_html_(self):
        if self.status == OPTIMAL:
            with option_context('display.max_rows', 10):
                html = ('<strong><em>Optimal</em> solution with objective '
                        'value {:.3f}</strong><br>{}'
                        .format(self.objective_value,
                                self.to_frame()._repr_html_()))
        else:
            html = '<strong><em>{}</em> solution</strong>'.format(self.status)
        return html

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
        warn("let Model.optimize create a solution instance,"
             " don't update yourself", DeprecationWarning)
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

    def to_frame(self):
        """Return the fluxes and reduced costs as a data frame"""
        return DataFrame({'fluxes': self.fluxes,
                          'reduced_costs': self.reduced_costs})


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
        Initialize a `LegacySolution` from an objective value.

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
            return "<LegacySolution {0:s} at 0x{1:x}>".format(
                self.status, id(self))
        return "<LegacySolution {0:.3f} at 0x{1:x}>".format(
            self.f, id(self))

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


def get_solution(model, reactions=None, metabolites=None, raise_error=False):
    """
    Generate a solution representation of the current solver state.

    Parameters
    ---------
    model : cobra.Model
        The model whose reactions to retrieve values for.
    reactions : list, optional
        An iterable of `cobra.Reaction` objects. Uses `model.reactions` by
        default.
    metabolites : list, optional
        An iterable of `cobra.Metabolite` objects. Uses `model.metabolites` by
        default.
    raise_error : bool
        If true, raise an OptimizationError if solver status is not optimal.

    Returns
    -------
    cobra.Solution

    Note
    ----
    This is only intended for the `optlang` solver interfaces and not the
    legacy solvers.
    """
    check_solver_status(model.solver.status, raise_error=raise_error)
    if reactions is None:
        reactions = model.reactions
    if metabolites is None:
        metabolites = model.metabolites

    rxn_index = list()
    fluxes = empty(len(reactions))
    reduced = empty(len(reactions))
    var_primals = model.solver.primal_values
    shadow = empty(len(metabolites))
    if model.solver.is_integer:
        reduced.fill(nan)
        shadow.fill(nan)
        for (i, rxn) in enumerate(reactions):
            rxn_index.append(rxn.id)
            fluxes[i] = var_primals[rxn.id] - var_primals[rxn.reverse_id]
        met_index = [met.id for met in metabolites]
    else:
        var_duals = model.solver.reduced_costs
        for (i, rxn) in enumerate(reactions):
            forward = rxn.id
            reverse = rxn.reverse_id
            rxn_index.append(forward)
            fluxes[i] = var_primals[forward] - var_primals[reverse]
            reduced[i] = var_duals[forward] - var_duals[reverse]
        met_index = list()
        constr_duals = model.solver.shadow_prices
        for (i, met) in enumerate(metabolites):
            met_index.append(met.id)
            shadow[i] = constr_duals[met.id]
    return Solution(model.solver.objective.value, model.solver.status,
                    Series(index=rxn_index, data=fluxes, name="fluxes"),
                    Series(index=rxn_index, data=reduced,
                           name="reduced_costs"),
                    Series(index=met_index, data=shadow, name="shadow_prices"))
