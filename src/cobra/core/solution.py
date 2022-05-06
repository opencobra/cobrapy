# -*- coding: utf-8 -*-

"""Provide unified interfaces to optimization solutions."""

from __future__ import absolute_import

import logging
from builtins import object, super

from numpy import empty, nan
from optlang.interface import OPTIMAL
from pandas import DataFrame, Series, option_context

from cobra.util.solver import check_solver_status


__all__ = ("Solution", "get_solution")

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
    """

    def __init__(
        self,
        objective_value,
        status,
        fluxes,
        reduced_costs=None,
        shadow_prices=None,
        **kwargs
    ):
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
        return "<Solution {0:.3f} at 0x{1:x}>".format(self.objective_value, id(self))

    def _repr_html_(self):
        if self.status == OPTIMAL:
            with option_context("display.max_rows", 10):
                html = (
                    "<strong><em>Optimal</em> solution with objective "
                    "value {:.3f}</strong><br>{}".format(
                        self.objective_value, self.to_frame()._repr_html_()
                    )
                )
        else:
            html = "<strong><em>{}</em> solution</strong>".format(self.status)
        return html

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

    def to_frame(self):
        """Return the fluxes and reduced costs as a data frame"""
        return DataFrame({"fluxes": self.fluxes, "reduced_costs": self.reduced_costs})


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
    return Solution(
        model.solver.objective.value,
        model.solver.status,
        Series(index=rxn_index, data=fluxes, name="fluxes"),
        Series(index=rxn_index, data=reduced, name="reduced_costs"),
        Series(index=met_index, data=shadow, name="shadow_prices"),
    )
