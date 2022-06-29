"""Provide functions for phenotype phase plane analysis."""

from itertools import product
from typing import TYPE_CHECKING, Dict, List, Optional, Union

import numpy as np
import pandas as pd
from optlang.interface import OPTIMAL

from ..exceptions import OptimizationError
from ..util import solver as sutil
from .helpers import normalize_cutoff
from .variability import flux_variability_analysis as fva


if TYPE_CHECKING:
    from optlang.interface import Objective

    from cobra import Model, Reaction


def production_envelope(
    model: "Model",
    reactions: List["Reaction"],
    objective: Union[Dict, "Objective", None] = None,
    carbon_sources: Optional[List["Reaction"]] = None,
    points: int = 20,
    threshold: Optional[float] = None,
) -> pd.DataFrame:
    """Calculate the objective value conditioned on all flux combinations.

    The production envelope can be used to analyze a model's ability to
    produce a given compound conditional on the fluxes for another set of
    reactions, such as the uptake rates. The model is alternately optimized
    with respect to minimizing and maximizing the objective and the
    obtained fluxes are recorded. Ranges to compute production is set to the
    effective bounds, i.e., the minimum / maximum fluxes that can be
    obtained given current reaction bounds.

    Parameters
    ----------
    model : cobra.Model
        The model to compute the production envelope for.
    reactions : list of cobra.Reaction
        A list of reaction objects.
    objective : dict or cobra.Model.objective, optional
        The objective (reaction) to use for the production envelope. Use the
        model's current objective if left missing (default None).
    carbon_sources : list of cobra.Reaction, optional
       One or more reactions that are the source of carbon for computing
       carbon (mol carbon in output over mol carbon in input) and mass
       yield (gram product over gram output). Only objectives with a carbon
       containing input and output metabolite is supported. Will identify
       active carbon sources in the medium if none are specified
       (default None).
    points : int, optional
       The number of points to calculate production for (default 20).
    threshold : float, optional
        A cut-off under which flux values will be considered to be zero.
        If not specified, it defaults to `model.tolerance` (default None).

    Returns
    -------
    pandas.DataFrame
        A DataFrame with fixed columns as:
        - carbon_source        : identifiers of carbon exchange reactions
        - flux_maximum         : maximum objective flux
        - flux_minimum         : minimum objective flux
        - carbon_yield_maximum : maximum yield of a carbon source
        - carbon_yield_minimum : minimum yield of a carbon source
        - mass_yield_maximum   : maximum mass yield of a carbon source
        - mass_yield_minimum   : minimum mass yield of a carbon source

        and variable columns (for each input `reactions`) as:
        - reaction_id          : flux at each given point

    Raises
    ------
    ValueError
        If model's objective is comprised of multiple reactions.

    Examples
    --------
    >>> import cobra.io
    >>> from cobra.flux_analysis import production_envelope
    >>> model = cobra.io.load_model("textbook")
    >>> production_envelope(model, ["EX_glc__D_e", "EX_o2_e"])
        carbon_source  flux_minimum  carbon_yield_minimum  mass_yield_minimum ...
    0     EX_glc__D_e           0.0                   0.0                 NaN ...
    1     EX_glc__D_e           0.0                   0.0                 NaN ...
    2     EX_glc__D_e           0.0                   0.0                 NaN ...
    3     EX_glc__D_e           0.0                   0.0                 NaN ...
    4     EX_glc__D_e           0.0                   0.0                 NaN ...
    ..            ...           ...                   ...                 ... ...
    395   EX_glc__D_e           NaN                   NaN                 NaN ...
    396   EX_glc__D_e           NaN                   NaN                 NaN ...
    397   EX_glc__D_e           NaN                   NaN                 NaN ...
    398   EX_glc__D_e           NaN                   NaN                 NaN ...
    399   EX_glc__D_e           NaN                   NaN                 NaN ...

    [400 rows x 9 columns]

    """
    reactions = model.reactions.get_by_any(reactions)
    objective = model.solver.objective if objective is None else objective
    data = {}

    if carbon_sources is None:
        c_input = _find_carbon_sources(model)
    else:
        c_input = model.reactions.get_by_any(carbon_sources)

    if c_input is None:
        data["carbon_source"] = None
    elif hasattr(c_input, "id"):
        data["carbon_source"] = c_input.id
    else:
        data["carbon_source"] = ", ".join(rxn.id for rxn in c_input)

    threshold = normalize_cutoff(model, threshold)

    size = points ** len(reactions)

    for direction in ("minimum", "maximum"):
        data[f"flux_{direction}"] = np.full(size, np.nan, dtype=float)
        data[f"carbon_yield_{direction}"] = np.full(size, np.nan, dtype=float)
        data[f"mass_yield_{direction}"] = np.full(size, np.nan, dtype=float)

    grid = pd.DataFrame(data)

    with model:
        model.objective = objective
        objective_reactions = list(sutil.linear_reaction_coefficients(model))

        if len(objective_reactions) != 1:
            raise ValueError(
                "Cannot calculate yields for objectives with multiple reactions."
            )
        c_output = objective_reactions[0]
        min_max = fva(model, reactions, fraction_of_optimum=0)
        min_max[min_max.abs() < threshold] = 0.0
        points = list(
            product(
                *[
                    np.linspace(
                        min_max.at[rxn.id, "minimum"],
                        min_max.at[rxn.id, "maximum"],
                        points,
                        endpoint=True,
                    )
                    for rxn in reactions
                ]
            )
        )
        tmp = pd.DataFrame(points, columns=[rxn.id for rxn in reactions])
        grid = pd.concat([grid, tmp], axis=1, copy=False)
        _add_envelope(model, reactions, grid, c_input, c_output, threshold)

    return grid


def _add_envelope(
    model: "Model",
    reactions: List["Reaction"],
    grid: pd.DataFrame,
    c_input: List["Reaction"],
    c_output: List["Reaction"],
    threshold: float,
) -> None:
    """Add a production envelope based on the parameters provided.

    Parameters
    ----------
    model : cobra.Model
        The model to operate on.
    reactions : list of cobra.Reaction
        The input reaction objects.
    grid : pandas.DataFrame
        The DataFrame containing all the data regarding the operation.
    c_input : list of cobra.Reaction
        The list of reaction objects acting as carbon inputs.
    c_output : list of cobra.Reaction
        The list of reaction objects acting as carbon outputs.

    """
    if c_input is not None:
        input_components = [_reaction_elements(rxn) for rxn in c_input]
        output_components = _reaction_elements(c_output)

        try:
            input_weights = [_reaction_weight(rxn) for rxn in c_input]
            output_weight = _reaction_weight(c_output)
        except ValueError:
            input_weights = []
            output_weight = []
    else:
        input_components = []
        output_components = []
        input_weights = []
        output_weight = []

    for direction in ("minimum", "maximum"):
        with model:
            model.objective_direction = direction

            for i in range(len(grid)):
                with model:
                    for rxn in reactions:
                        point = grid.at[i, rxn.id]
                        rxn.bounds = point, point
                    obj_val = model.slim_optimize()

                    if model.solver.status != OPTIMAL:
                        continue

                    grid.at[i, f"flux_{direction}"] = (
                        0.0 if np.abs(obj_val) < threshold else obj_val
                    )

                    if c_input is not None:
                        grid.at[i, f"carbon_yield_{direction}"] = _total_yield(
                            [rxn.flux for rxn in c_input],
                            input_components,
                            obj_val,
                            output_components,
                        )
                        grid.at[i, f"mass_yield_{direction}"] = _total_yield(
                            [rxn.flux for rxn in c_input],
                            input_weights,
                            obj_val,
                            output_weight,
                        )


def _total_yield(
    input_fluxes: List[float],
    input_elements: List[float],
    output_flux: List[float],
    output_elements: List[float],
) -> float:
    """Compute total output per input unit.

    Units are typically mol carbon atoms or gram of source and product.

    Parameters
    ----------
    input_fluxes : list of float
        A list of input reaction fluxes in the same order as the
        `input_components`.
    input_elements : list of float
        A list of reaction components which are in turn list of numbers.
    output_flux : float
        The output flux value.
    output_elements : list
        A list of stoichiometrically weighted output reaction components.

    Returns
    -------
    float
        The ratio between output (mol carbon atoms or grams of product) and
        input (mol carbon atoms or grams of source compounds). If input flux
        of carbon sources is zero then numpy.nan is returned.

    """
    carbon_input_flux = sum(
        _total_components_flux(flux, components, consumption=True)
        for flux, components in zip(input_fluxes, input_elements)
    )
    carbon_output_flux = _total_components_flux(
        output_flux, output_elements, consumption=False
    )
    try:
        return carbon_output_flux / carbon_input_flux
    except ZeroDivisionError:
        return np.nan


def _reaction_elements(reaction: "Reaction") -> List[float]:
    """Split metabolites into atoms times their stoichiometric coefficients.

    Parameters
    ----------
    reaction : cobra.Reaction
        The reaction whose metabolite components are desired.

    Returns
    -------
    list of float
        Each of the reaction's metabolites' desired carbon elements (if any)
        times that metabolite's stoichiometric coefficient.

    """
    c_elements = [
        coeff * met.elements.get("C", 0) for met, coeff in reaction.metabolites.items()
    ]
    return [elem for elem in c_elements if elem != 0]


def _reaction_weight(reaction: "Reaction") -> List[float]:
    """Return the metabolite weight times its stoichiometric coefficient.

    Parameters
    ----------
    reaction : cobra.Reaction
        The reaction whose metabolite component weights is desired.

    Returns
    -------
    list of float
        Each of reaction's metabolite components' weights.

    Raises
    ------
    ValueError
        If more than one metabolite comprises the `reaction`.

    """
    if len(reaction.metabolites) != 1:
        raise ValueError(
            "Reaction weight is only defined for single "
            "metabolite products or educts."
        )

    met, coeff = next(iter(reaction.metabolites.items()))
    return [coeff * met.formula_weight]


def _total_components_flux(
    flux: float, components: List[float], consumption: bool = True
) -> float:
    """Compute the total components consumption or production flux.

    Parameters
    ----------
    flux : float
        The reaction flux for the components.
    components : list of float
        List of stoichiometrically weighted components.
    consumption : bool, optional
        Whether to sum up consumption or production fluxes (default True).

    Returns
    -------
    float
        The total consumption or production flux of `components`.

    """
    direction = 1 if consumption else -1
    c_flux = [elem * flux * direction for elem in components]
    return sum([flux for flux in c_flux if flux > 0])


def _find_carbon_sources(model: "Model") -> List["Reaction"]:
    """Find all active carbon source reactions.

    Parameters
    ----------
    model : Model
        The model whose active carbon sources need to found.

    Returns
    -------
    list of cobra.Reaction
       The medium reactions with carbon input flux.

    """
    try:
        model.slim_optimize(error_value=None)
    except OptimizationError:
        return []

    reactions = model.reactions.get_by_any(list(model.medium))
    reactions_fluxes = [
        (
            rxn,
            _total_components_flux(rxn.flux, _reaction_elements(rxn), consumption=True),
        )
        for rxn in reactions
    ]
    return [rxn for rxn, c_flux in reactions_fluxes if c_flux > 0]
