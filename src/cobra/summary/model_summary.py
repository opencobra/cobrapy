"""Provide the model summary class."""


import logging
from operator import attrgetter
from typing import TYPE_CHECKING, Dict, List, Optional, Union

import pandas as pd

from cobra.core import Reaction
from cobra.flux_analysis import flux_variability_analysis, pfba
from cobra.summary import Summary
from cobra.util.solver import linear_reaction_coefficients


if TYPE_CHECKING:
    from cobra.core import Metabolite, Model, Solution


logger = logging.getLogger(__name__)


class ModelSummary(Summary):
    """
    Define the model summary.

    Attributes
    ----------
    uptake_flux : pandas.DataFrame
        A pandas DataFrame of only the uptake fluxes.
    secretion_flux : pandas.DataFrame
        A pandas DataFrame of only the consuming fluxes.
    See Also
    --------
    Summary : Parent that defines further attributes.
    MetaboliteSummary
    ReactionSummary

    """

    def __init__(
        self,
        *,
        model: "Model",
        solution: Optional["Solution"] = None,
        fva: Optional[Union[float, pd.DataFrame]] = None,
        **kwargs,
    ):
        """
        Initialize a model summary.

        Parameters
        ----------
        model : cobra.Model
            The metabolic model for which to generate a metabolite summary.
        solution : cobra.Solution, optional
            A previous model solution to use for generating the summary. If
            ``None``, the summary method will generate a parsimonious flux
            distribution (default None).
        fva : pandas.DataFrame or float, optional
            Whether or not to include flux variability analysis in the output.
            If given, `fva` should either be a previous FVA solution matching the
            model or a float between 0 and 1 representing the fraction of the
            optimum objective to be searched (default None).

        Other Parameters
        ----------------
        kwargs :
            Further keyword arguments are passed on to the parent class.

        See Also
        --------
        Summary : Parent that has further default parameters.
        MetaboliteSummary
        ReactionSummary

        """
        super().__init__(**kwargs)
        self._objective = None
        self._objective_value = None
        self._boundary = None
        self._boundary_metabolites = None
        self.uptake_flux: Optional[pd.DataFrame] = None
        self.secretion_flux: Optional[pd.DataFrame] = None
        self._generate(model, solution, fva)

    def _generate(
        self,
        model: "Model",
        solution: Optional["Solution"],
        fva: Optional[Union[float, pd.DataFrame]],
    ) -> None:
        """
        Prepare the data for the summary instance.

        Parameters
        ----------
        model : cobra.Model
            The metabolic model for which to generate a metabolite summary.
        solution : cobra.Solution, optional
            A previous model solution to use for generating the summary. If
            ``None``, the summary method will generate a parsimonious flux
            distribution.
        fva : pandas.DataFrame or float, optional
            Whether or not to include flux variability analysis in the output.
            If given, `fva` should either be a previous FVA solution matching the
            model or a float between 0 and 1 representing the fraction of the
            optimum objective to be searched.

        """
        super()._generate(model=model, solution=solution, fva=fva)

        coefficients = linear_reaction_coefficients(model)
        if solution is None:
            logger.info("Generating new parsimonious flux distribution.")
            solution = pfba(model)

        if isinstance(fva, float):
            logger.info("Performing flux variability analysis.")
            fva = flux_variability_analysis(
                model=model,
                reaction_list=model.boundary,
                fraction_of_optimum=fva,
            )
        if coefficients:
            self._objective: Dict["Reaction", float] = {
                rxn.copy(): coef for rxn, coef in coefficients.items()
            }
            self._objective_value: float = sum(
                solution[rxn.id] * coef for rxn, coef in self._objective.items()
            )
        else:
            logger.warning(
                "Non-linear or non-reaction model objective. Falling back to minimal "
                "display."
            )
            self._objective = {
                Reaction(id="Expression", name="Expression"): float("nan")
            }
            self._objective_value: float = float("nan")
        self._boundary: List["Reaction"] = [
            rxn.copy() for rxn in sorted(model.boundary, key=attrgetter("id"))
        ]
        self._boundary_metabolites: List["Metabolite"] = [
            met.copy() for rxn in self._boundary for met in rxn.metabolites
        ]
        flux = pd.DataFrame(
            data=[
                (rxn.id, met.id, rxn.get_coefficient(met.id), solution[rxn.id])
                for rxn, met in zip(self._boundary, self._boundary_metabolites)
            ],
            columns=["reaction", "metabolite", "factor", "flux"],
            index=[r.id for r in self._boundary],
        )
        # Scale fluxes by stoichiometric coefficient.
        flux["flux"] *= flux["factor"]

        if fva is not None:
            flux = flux.join(fva)
            view = flux[["flux", "minimum", "maximum"]]
            # Set fluxes below model tolerance to zero.
            flux[["flux", "minimum", "maximum"]] = view.where(
                view.abs() >= model.tolerance, 0
            )
            # Create the scaled compound flux.
            flux[["minimum", "maximum"]] = flux[["minimum", "maximum"]].mul(
                flux["factor"], axis=0
            )
            # Negative factors invert the minimum/maximum relationship.
            negative = flux["factor"] < 0
            tmp = flux.loc[negative, "maximum"]
            flux.loc[negative, "maximum"] = flux.loc[negative, "minimum"]
            flux.loc[negative, "minimum"] = tmp
            # Add zero to turn negative zero into positive zero for nicer display later.
            flux[["flux", "minimum", "maximum"]] += 0
        else:
            # Set fluxes below model tolerance to zero.
            flux.loc[flux["flux"].abs() < model.tolerance, "flux"] = 0
            # Add zero to turn negative zero into positive zero for nicer display later.
            flux["flux"] += 0

        # Create production table from producing fluxes or zero fluxes where the
        # metabolite is a product in the reaction.
        is_produced = (flux["flux"] > 0) | ((flux["flux"] == 0) & (flux["factor"] > 0))
        if fva is not None:
            self.uptake_flux = flux.loc[
                is_produced, ["flux", "minimum", "maximum", "reaction", "metabolite"]
            ].copy()
        else:
            self.uptake_flux = flux.loc[
                is_produced, ["flux", "reaction", "metabolite"]
            ].copy()

        # Create consumption table from consuming fluxes or zero fluxes where the
        # metabolite is a substrate in the reaction.
        is_consumed = (flux["flux"] < 0) | ((flux["flux"] == 0) & (flux["factor"] < 0))
        if fva is not None:
            self.secretion_flux = flux.loc[
                is_consumed, ["flux", "minimum", "maximum", "reaction", "metabolite"]
            ].copy()
        else:
            self.secretion_flux = flux.loc[
                is_consumed, ["flux", "reaction", "metabolite"]
            ].copy()

        self._flux = flux

    def _display_flux(
        self, frame: pd.DataFrame, names: bool, element: str, threshold: float
    ) -> pd.DataFrame:
        """
        Transform a flux data frame for display.

        Parameters
        ----------
        frame : pandas.DataFrame
            Either the producing or the consuming fluxes.
        names : bool
            Whether or not elements should be displayed by their common names.
        element : str
            The atomic element to summarize fluxes for.
        threshold : float
            Hide fluxes below the threshold from being displayed.

        Returns
        -------
        pandas.DataFrame
            The transformed pandas DataFrame with flux percentages and reaction
            definitions.

        """
        if "minimum" in frame.columns and "maximum" in frame.columns:
            frame = frame.loc[
                (frame["flux"].abs() >= threshold)
                | (frame["minimum"].abs() >= threshold)
                | (frame["maximum"].abs() >= threshold),
                :,
            ].copy()
        else:
            frame = frame.loc[frame["flux"].abs() >= threshold, :].copy()

        metabolites = {m.id: m for m in self._boundary_metabolites}

        element_num = f"{element}-Number"
        frame[element_num] = [
            metabolites[met_id].elements.get(element, 0)
            for met_id in frame["metabolite"]
        ]
        element_percent = f"{element}-Flux"
        frame[element_percent] = frame[element_num] * frame["flux"].abs()
        total = frame[element_percent].sum()
        if total > 0.0:
            frame[element_percent] /= total
        frame[element_percent] = [f"{x:.2%}" for x in frame[element_percent]]

        if names:
            frame["metabolite"] = [
                metabolites[met_id].name for met_id in frame["metabolite"]
            ]

        if "minimum" in frame.columns and "maximum" in frame.columns:
            frame["range"] = list(
                frame[["minimum", "maximum"]].itertuples(index=False, name=None)
            )
            return frame[
                [
                    "metabolite",
                    "reaction",
                    "flux",
                    "range",
                    element_num,
                    element_percent,
                ]
            ]
        else:
            return frame[
                ["metabolite", "reaction", "flux", element_num, element_percent]
            ]

    @staticmethod
    def _string_table(frame: pd.DataFrame, float_format: str, column_width: int) -> str:
        """
        Create a pretty string representation of the data frame.

        Parameters
        ----------
        frame : pandas.DataFrame
            A pandas DataFrame of fluxes.
        float_format : str
            Format string for floats.
        column_width : int
            The maximum column width for each row.

        Returns
        -------
        str
            The data frame formatted as a pretty string.

        """
        frame.columns = [header.title() for header in frame.columns]
        return frame.to_string(
            header=True,
            index=False,
            na_rep="",
            formatters={
                "Flux": f"{{:{float_format}}}".format,
                "Range": lambda pair: f"[{pair[0]:{float_format}}; "
                f"{pair[1]:{float_format}}]",
            },
            max_colwidth=column_width,
        )

    @staticmethod
    def _html_table(frame: pd.DataFrame, float_format: str) -> str:
        """
        Create an HTML representation of the data frame.

        Parameters
        ----------
        frame : pandas.DataFrame
            A pandas DataFrame of fluxes.
        float_format : str
            Format string for floats.

        Returns
        -------
        str
            The data frame formatted as HTML.

        """
        frame.columns = [header.title() for header in frame.columns]
        return frame.to_html(
            header=True,
            index=False,
            na_rep="",
            formatters={
                "Flux": f"{{:{float_format}}}".format,
                "Range": lambda pair: f"[{pair[0]:{float_format}}; "
                f" {pair[1]:{float_format}}]",
            },
        )

    def _string_objective(self, names: bool) -> str:
        """
        Return a string representation of the objective.

        Parameters
        ----------
        names : bool, optional
            Whether or not elements should be displayed by their common names.

        Returns
        -------
        str
            The objective expression and value as a string.

        """
        if names:
            objective = " + ".join(
                [f"{coef} {rxn.name}" for rxn, coef in self._objective.items()]
            )
        else:
            objective = " + ".join(
                [f"{coef} {rxn.id}" for rxn, coef in self._objective.items()]
            )
        return f"{objective} = {self._objective_value}"

    def to_string(
        self,
        names: bool = False,
        element: str = "C",
        threshold: Optional[float] = None,
        float_format: str = ".4G",
        column_width: int = 79,
    ) -> str:
        """
        Return a pretty string representation of the model summary.

        Parameters
        ----------
        names : bool, optional
            Whether or not elements should be displayed by their common names
            (default False).
        element : str, optional
            The atomic element to summarize uptake and secretion for (default 'C').
        threshold : float, optional
            Hide fluxes below the threshold from being displayed. If no value is
            given, the model tolerance is used (default None).
        float_format : str, optional
            Format string for floats (default '.4G').
        column_width : int, optional
            The maximum column width for each row (default 79).

        Returns
        -------
        str
            The summary formatted as a pretty string.

        """
        threshold = self._normalize_threshold(threshold)

        objective = self._string_objective(names)

        uptake = self._string_table(
            self._display_flux(self.uptake_flux, names, element, threshold),
            float_format,
            column_width,
        )

        secretion = self._string_table(
            self._display_flux(self.secretion_flux, names, element, threshold),
            float_format,
            column_width,
        )

        return (
            f"Objective\n"
            f"=========\n"
            f"{objective}\n\n"
            f"Uptake\n"
            f"------\n"
            f"{uptake}\n\n"
            f"Secretion\n"
            f"---------\n"
            f"{secretion}\n"
        )

    def to_html(
        self,
        names: bool = False,
        element: str = "C",
        threshold: Optional[float] = None,
        float_format: str = ".4G",
    ) -> str:
        """
        Return a rich HTML representation of the model summary.

        Parameters
        ----------
        names : bool, optional
            Whether or not elements should be displayed by their common names
            (default False).
        element : str, optional
            The atomic element to summarize uptake and secretion for (default 'C').
        threshold : float, optional
            Hide fluxes below the threshold from being displayed. If no value is
            given, the model tolerance is used (default None).
        float_format : str, optional
            Format string for floats (default '.4G').

        Returns
        -------
        str
            The summary formatted as HTML.

        """
        threshold = self._normalize_threshold(threshold)

        objective = self._string_objective(names)

        uptake = self._html_table(
            self._display_flux(self.uptake_flux, names, element, threshold),
            float_format,
        )

        secretion = self._html_table(
            self._display_flux(self.secretion_flux, names, element, threshold),
            float_format,
        )

        return (
            f"<h3>Objective</h3>"
            f"<p>{objective}</p>"
            f"<h4>Uptake</h4>"
            f"{uptake}"
            f"<h4>Secretion</h4>"
            f"{secretion}"
        )
