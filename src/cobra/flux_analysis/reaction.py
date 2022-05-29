"""Provide functions for analyzing / creating objective functions."""

from operator import attrgetter
from typing import TYPE_CHECKING, Dict, Union
from warnings import warn

from ..core import Reaction


if TYPE_CHECKING:
    from cobra import Model


def assess(
    model: "Model",
    reaction: Union[Reaction, str],
    flux_coefficient_cutoff: float = 0.001,
    solver: str = None,
) -> Union[bool, Dict]:
    """Assess production capacity.

    Assess the capacity of the `model` to produce the precursors for the
    `reaction` and absorb the production of the `reaction` while the
    `reaction` is operating at, or above, the specified
    `flux_coefficient_cutoff`.

    .. deprecated:: 0.10.0a1
              `solver` argument will be removed in cobrapy 1.0.0, it is
              replaced by globally setting the solver via
              cobra.Configuration .

    Parameters
    ----------
    model : cobra.Model
        The cobra model to assess production capacity for.
    reaction : str or cobra.Reaction
        The reaction to assess.
    flux_coefficient_cutoff : float, optional
        The minimum flux that reaction must carry to be considered active
        (default 0.001).
    solver : str, optional
        The solver name. If None, the default solver will be used
        (default None).

    Returns
    -------
    bool or dict
        True if the model can produce the precursors and absorb the products
        for the reaction operating at, or above, `flux_coefficient_cutoff`.
        Otherwise, a dictionary of {'precursor': Status, 'product': Status},
        where 'Status' is the results from `assess_precursors` and
        `assess_products`, respectively.

    """
    reaction = model.reactions.get_by_any(reaction)[0]
    with model as m:
        m.objective = reaction
        if m.slim_optimize(error_value=0.0) >= flux_coefficient_cutoff:
            return True
        else:
            results = {}
            results["precursors"] = assess_component(
                model, reaction, "reactants", flux_coefficient_cutoff
            )
            results["products"] = assess_component(
                model, reaction, "products", flux_coefficient_cutoff
            )
            return results


def assess_component(
    model: "Model",
    reaction: Union[Reaction, str],
    side: str,
    flux_coefficient_cutoff: float = 0.001,
    solver: str = None,
) -> Union[bool, Dict]:
    """Assess production capacity of components.

    Assess the ability of the `model` to provide sufficient precursors,
    or absorb products, for a `reaction` operating at, or beyond,
    the specified `flux_coefficient_cutoff`.

    .. deprecated:: 0.10.0a1
              `solver` argument will be removed in cobrapy 1.0.0, it is
              replaced by globally setting the solver via
              cobra.Configuration .

    Parameters
    ----------
    model : cobra.Model
        The cobra model to assess production capacity for.
    reaction : str or cobra.Reaction
        The reaction to assess.
    side : {"products", "reactants"}
        The side of the reaction to assess.
    flux_coefficient_cutoff : float, optional
        The minimum flux that reaction must carry to be considered active
        (default 0.001).
    solver : str, optional
        The solver name. If None, the default solver will be used
        (default None).

    Returns
    -------
    bool or dict
        True if the precursors can be simultaneously produced at the
        specified cutoff. False, if the model has the capacity to produce
        each individual precursor at the specified threshold  but not all
        precursors at the required level simultaneously. Otherwise, a
        dictionary of the required and the produced fluxes for each
        reactant that is not produced in sufficient quantities.

    """
    reaction = model.reactions.get_by_any(reaction)[0]
    result_key = {"reactants": "produced", "products": "capacity"}[side]
    get_components = attrgetter(side)
    with model as m:
        m.objective = reaction
        if m.slim_optimize(error_value=0.0) >= flux_coefficient_cutoff:
            return True
        simulation_results = {}
        # build the demand reactions and add all at once
        demand_reactions = {}
        for component in get_components(reaction):
            coeff = reaction.metabolites[component]
            demand = m.add_boundary(component, type="demand")
            demand.metabolites[component] = coeff
            demand_reactions[demand] = (component, coeff)
        # First assess whether all precursors can be produced simultaneously
        joint_demand = Reaction("joint_demand")
        for demand_reaction in demand_reactions:
            joint_demand += demand_reaction
        m.add_reactions([joint_demand])
        m.objective = joint_demand
        if m.slim_optimize(error_value=0.0) >= flux_coefficient_cutoff:
            return True

        # Otherwise assess the ability of the model to produce each precursor
        # individually.  Now assess the ability of the model to produce each
        # reactant for a reaction
        for demand_reaction, (component, coeff) in demand_reactions.items():
            # Calculate the maximum amount of the
            with m:
                m.objective = demand_reaction
                flux = m.slim_optimize(error_value=0.0)
            # metabolite that can be produced.
            if flux_coefficient_cutoff > flux:
                # Scale the results to a single unit
                simulation_results.update(
                    {
                        component: {
                            "required": flux_coefficient_cutoff / abs(coeff),
                            result_key: flux / abs(coeff),
                        }
                    }
                )
        if len(simulation_results) == 0:
            simulation_results = False
        return simulation_results


def _optimize_or_value(model: "Model", value: float = 0.0) -> float:
    """Perform quick optimization and return the objective value.

    The objective value is returned at `value` error value.

    .. deprecated:: 0.22.0
              `_optimize_or_value` will be removed in cobrapy 1.0.0, its
              functionality can be achieved by using
              `cobra.Model.slim_optimize`.

    Parameters
    ----------
    model: cobra.Model
        The model to optimize.
    value: float
        The error value to consider.

    Returns
    -------
    float
        The optimized value of the model's objective.

    """
    return model.slim_optimize(error_value=value)


def assess_precursors(
    model: "Model",
    reaction: Reaction,
    flux_coefficient_cutoff: float = 0.001,
    solver: str = None,
) -> Union[bool, Dict]:
    """Assess production capacity of precursors.

    Assess the ability of the model to provide sufficient precursors for
    a reaction operating at, or beyond, the `flux_coefficient_cutoff`.

    .. deprecated:: 0.7.0
              `assess_precursors` will be removed in cobrapy 1.0.0, it is
              replaced by `assess_component`.

    Parameters
    ----------
    model : cobra.Model
        The cobra model to assess production capacity for.
    reaction : str or cobra.Reaction
        The reaction to assess.
    flux_coefficient_cutoff : float, optional
        The minimum flux that reaction must carry to be considered active
        (default 0.001).
    solver : str, optional
        The solver name. If None, the default solver will be used
        (default None).

    Returns
    -------
    bool or dict
        True if the precursors can be simultaneously produced at the
        specified cutoff. False, if the model has the capacity to produce
        each individual precursor at the specified threshold but not all
        precursors at the required level simultaneously. Otherwise, a
        dictionary of the required and the produced fluxes for each
        reactant that is not produced in sufficient quantities.

    """
    warn("Use cobra.sampling.reaction.assess_component() instead.", DeprecationWarning)
    return assess_component(
        model, reaction, "reactants", flux_coefficient_cutoff, solver
    )


def assess_products(
    model: "Model",
    reaction: Reaction,
    flux_coefficient_cutoff: float = 0.001,
    solver: str = None,
) -> Union[bool, Dict]:
    """Assesses absorption capacity of products.

    Assess the ability of the model to to absorb the products of a reaction
    at a given flux rate. Useful for identifying which components might be
    blocking a reaction from achieving a specific flux rate.

    .. deprecated:: 0.7.0
              `assess_products` will be removed in cobrapy 1.0.0, it is
              replaced by `assess_component`.

    Parameters
    ----------
    model : cobra.Model
        The cobra model to assess production capacity for.
    reaction : str or cobra.Reaction
        The reaction to assess.
    flux_coefficient_cutoff :  float, optional
        The minimum flux that reaction must carry to be considered active
        (default 0.001).
    solver : str, optional
        The solver name. If None, the default solver will be used
        (default None).

    Returns
    -------
    bool or dict
        True if the model has the capacity to absorb all the reaction
        products being simultaneously given the specified cutoff. False,
        if the model has the capacity to absorb each individual product but
        not all products at the required level simultaneously. Otherwise, a
        dictionary of the required and the capacity fluxes for each product
        that is not absorbed in sufficient quantities.

    """
    warn("Use cobra.sampling.reaction.assess_component() instead.", DeprecationWarning)
    return assess_component(
        model, reaction, "products", flux_coefficient_cutoff, solver
    )
