"""Provide the metabolite summary class."""


import logging
from operator import attrgetter
from textwrap import shorten
from typing import TYPE_CHECKING, List, Optional, Union

import pandas as pd

from cobra.flux_analysis import flux_variability_analysis, pfba
from cobra.summary import Summary


if TYPE_CHECKING:
    from cobra.core import Metabolite, Model, Reaction, Solution


logger = logging.getLogger(__name__)


class MetaboliteSummary(Summary):
    """
    Define the metabolite summary.

    Attributes
    ----------
    producing_flux : pandas.DataFrame
        A pandas DataFrame of only the producing fluxes.
    consuming_flux : pandas.DataFrame
        A pandas DataFrame of only the consuming fluxes.

    See Also
    --------
    Summary : Parent that defines further attributes.
    ReactionSummary
    ModelSummary

    """

    def __init__(
        self,
        *,
        metabolite: "Metabolite",
        model: "Model",
        solution: Optional["Solution"] = None,
        fva: Optional[Union[float, pd.DataFrame]] = None,
        **kwargs,
    ) -> None:
        """
        Initialize a metabolite summary.

        Parameters
        ----------
        metabolite : cobra.Metabolite
            The metabolite object whose summary we intend to get.
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
        ReactionSummary
        ModelSummary

        """
        super().__init__(**kwargs)
        self._metabolite = metabolite.copy()
        self._reactions: List["Reaction"] = [
            r.copy() for r in sorted(metabolite.reactions, key=attrgetter("id"))
        ]
        self.producing_flux: Optional[pd.DataFrame] = None
        self.consuming_flux: Optional[pd.DataFrame] = None
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

        if solution is None:
            logger.info("Generating new parsimonious flux distribution.")
            solution = pfba(model)

        if isinstance(fva, float):
            logger.info("Performing flux variability analysis.")
            fva = flux_variability_analysis(
                model=model,
                reaction_list=[r.id for r in self._reactions],
                fraction_of_optimum=fva,
            )

        # Create the basic flux table.
        flux = pd.DataFrame(
            data=[
                (
                    r.id,
                    solution[r.id],
                    r.get_coefficient(self._metabolite.id),
                )
                for r in self._reactions
            ],
            columns=["reaction", "flux", "factor"],
            index=[r.id for r in self._reactions],
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
            self.producing_flux = flux.loc[
                is_produced, ["flux", "minimum", "maximum", "reaction"]
            ].copy()
        else:
            self.producing_flux = flux.loc[is_produced, ["flux", "reaction"]].copy()
        production = self.producing_flux["flux"].abs()
        self.producing_flux["percent"] = production / production.sum()

        # Create consumption table from consuming fluxes or zero fluxes where the
        # metabolite is a substrate in the reaction.
        is_consumed = (flux["flux"] < 0) | ((flux["flux"] == 0) & (flux["factor"] < 0))
        if fva is not None:
            self.consuming_flux = flux.loc[
                is_consumed, ["flux", "minimum", "maximum", "reaction"]
            ].copy()
        else:
            self.consuming_flux = flux.loc[is_consumed, ["flux", "reaction"]].copy()
        consumption = self.consuming_flux["flux"].abs()
        self.consuming_flux["percent"] = consumption / consumption.sum()

        self._flux = flux

    def _display_flux(
        self, frame: pd.DataFrame, names: bool, threshold: float
    ) -> pd.DataFrame:
        """
        Transform a flux data frame for display.

        Parameters
        ----------
        frame : pandas.DataFrame
            Either the producing or the consuming fluxes.
        names : bool
            Whether or not elements should be displayed by their common names.
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
        reactions = {r.id: r for r in self._reactions}
        frame["definition"] = [
            reactions[rxn_id].build_reaction_string(names)
            for rxn_id in frame["reaction"]
        ]
        if "minimum" in frame.columns and "maximum" in frame.columns:
            frame["range"] = list(
                frame[["minimum", "maximum"]].itertuples(index=False, name=None)
            )
            return frame[["percent", "flux", "range", "reaction", "definition"]]
        else:
            return frame[["percent", "flux", "reaction", "definition"]]

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
                "Percent": "{:.2%}".format,
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
                "Percent": "{:.2%}".format,
                "Flux": f"{{:{float_format}}}".format,
                "Range": lambda pair: f"[{pair[0]:{float_format}}; "
                f" {pair[1]:{float_format}}]",
            },
        )

    def to_string(
        self,
        names: bool = False,
        threshold: Optional[float] = None,
        float_format: str = ".4G",
        column_width: int = 79,
    ) -> str:
        """
        Return a pretty string representation of the metabolite summary.

        Parameters
        ----------
        names : bool, optional
            Whether or not elements should be displayed by their common names
            (default False).
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

        if names:
            metabolite = shorten(
                self._metabolite.name, width=column_width, placeholder="..."
            )
        else:
            metabolite = shorten(
                self._metabolite.id, width=column_width, placeholder="..."
            )

        production = self._string_table(
            self._display_flux(self.producing_flux, names, threshold),
            float_format,
            column_width,
        )

        consumption = self._string_table(
            self._display_flux(self.consuming_flux, names, threshold),
            float_format,
            column_width,
        )

        return (
            f"{metabolite}\n"
            f"{'=' * len(metabolite)}\n"
            f"Formula: {self._metabolite.formula}\n\n"
            f"Producing Reactions\n"
            f"-------------------\n"
            f"{production}\n\n"
            f"Consuming Reactions\n"
            f"-------------------\n"
            f"{consumption}"
        )

    def to_html(
        self,
        names: bool = False,
        threshold: Optional[float] = None,
        float_format: str = ".4G",
    ) -> str:
        """
        Return a rich HTML representation of the metabolite summary.

        Parameters
        ----------
        names : bool, optional
            Whether or not elements should be displayed by their common names
            (default False).
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

        if names:
            metabolite = self._metabolite.name
        else:
            metabolite = self._metabolite.id

        production = self._html_table(
            self._display_flux(self.producing_flux, names, threshold),
            float_format,
        )

        consumption = self._html_table(
            self._display_flux(self.consuming_flux, names, threshold),
            float_format,
        )

        return (
            f"<h3>{metabolite}</h3>"
            f"<p>{self._metabolite.formula}</p>"
            f"<h4>Producing Reactions</h4>"
            f"{production}"
            f"<h4>Consuming Reactions</h4>"
            f"{consumption}"
        )
