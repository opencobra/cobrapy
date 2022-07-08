"""Provide unified interfaces to optimization solutions."""

import logging
from typing import TYPE_CHECKING, Iterable, Optional

import numpy as np
import pandas as pd
from optlang.interface import OPTIMAL

from ..util.solver import check_solver_status


if TYPE_CHECKING:
    from cobra import Metabolite, Model, Reaction


__all__ = ("Solution", "get_solution")

logger = logging.getLogger(__name__)


class Solution:
    """
    A unified interface to a `cobra.Model` optimization solution.

    Parameters
    ----------
    objective_value : float
        The (optimal) value for the objective function.
    status : str
        The solver status related to the solution.
    fluxes : pandas.Series
        Contains the reaction fluxes (primal values of variables).
    reduced_costs : pandas.Series
        Contains reaction reduced costs (dual values of variables)
        (default None).
    shadow_prices : pandas.Series
        Contains metabolite shadow prices (dual values of constraints)
        (default None).

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

    Notes
    -----
    Solution is meant to be constructed by `get_solution` please look at that
    function to fully understand the `Solution` class.

    """

    def __init__(
        self,
        objective_value: float,
        status: str,
        fluxes: pd.Series,
        reduced_costs: Optional[pd.Series] = None,
        shadow_prices: Optional[pd.Series] = None,
        **kwargs,
    ) -> None:
        """
        Initialize a `Solution` from its components.

        Other Parameters
        ----------------
        kwargs :
            Further keyword arguments are passed on to the parent class.

        """
        super().__init__(**kwargs)
        self.objective_value = objective_value
        self.status = status
        self.fluxes = fluxes
        self.reduced_costs = reduced_costs
        self.shadow_prices = shadow_prices

    def __repr__(self) -> str:
        """Return a string representation of the solution instance."""
        if self.status != OPTIMAL:
            return f"<Solution {self.status} at {id(self):#x}>"
        return f"<Solution {self.objective_value:.3f} at {id(self):#x}>"

    def _repr_html_(self) -> str:
        """Return a rich HTML representation of the solution."""
        if self.status == OPTIMAL:
            with pd.option_context("display.max_rows", 10):
                html = (
                    "<strong><em>Optimal</em> solution with objective "
                    f"value {self.objective_value:.3f}</strong><br>"
                    f"{self.to_frame()._repr_html_()}"
                )
        else:
            html = f"<strong><em>{self.status}</em> solution</strong>"
        return html

    def __getitem__(self, reaction_id: str) -> float:
        """
        Return the flux of a reaction.

        Parameters
        ----------
        reaction_id : str
            A model reaction ID.

        Returns
        -------
        float
            The flux of the reaction with ID `reaction_id`.

        """
        return self.fluxes[reaction_id]

    get_primal_by_id = __getitem__

    def to_frame(self) -> pd.DataFrame:
        """Return the fluxes and reduced costs as a pandas DataFrame.

        Returns
        -------
        pandas.DataFrame
            The fluxes and reduced cost.

        """
        return pd.DataFrame(
            {"fluxes": self.fluxes, "reduced_costs": self.reduced_costs}
        )


def get_solution(
    model: "Model",
    reactions: Optional[Iterable["Reaction"]] = None,
    metabolites: Optional[Iterable["Metabolite"]] = None,
    raise_error: bool = False,
) -> Solution:
    """
    Generate a solution representation of the current solver state.

    Parameters
    ---------
    model : cobra.Model
        The model whose reactions to retrieve values for.
    reactions : list, optional
        An iterable of `cobra.Reaction` objects. Uses `model.reactions`
        if None (default None).
    metabolites : list, optional
        An iterable of `cobra.Metabolite` objects. Uses `model.metabolites`
        if None (default None).
    raise_error : bool
        If True, raise an OptimizationError if solver status is not optimal
        (default False).

    Returns
    -------
    cobra.Solution

    """
    check_solver_status(model.solver.status, raise_error=raise_error)
    if reactions is None:
        reactions = model.reactions
    if metabolites is None:
        metabolites = model.metabolites

    rxn_index = []
    fluxes = np.empty(len(reactions))
    reduced = np.empty(len(reactions))
    var_primals = model.solver.primal_values
    shadow = np.empty(len(metabolites))
    if model.solver.is_integer:
        reduced.fill(np.nan)
        shadow.fill(np.nan)
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
        met_index = []
        constr_duals = model.solver.shadow_prices
        for (i, met) in enumerate(metabolites):
            met_index.append(met.id)
            shadow[i] = constr_duals[met.id]
    return Solution(
        objective_value=model.solver.objective.value,
        status=model.solver.status,
        fluxes=pd.Series(index=rxn_index, data=fluxes, name="fluxes"),
        reduced_costs=pd.Series(index=rxn_index, data=reduced, name="reduced_costs"),
        shadow_prices=pd.Series(index=met_index, data=shadow, name="shadow_prices"),
    )
