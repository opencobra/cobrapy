"""Provide the reaction summary class."""


import logging
from textwrap import dedent, shorten
from typing import TYPE_CHECKING, Optional, Union

import pandas as pd

from cobra.flux_analysis import flux_variability_analysis, pfba
from cobra.summary import Summary


if TYPE_CHECKING:
    from cobra import Model, Reaction, Solution


logger = logging.getLogger(__name__)


class ReactionSummary(Summary):
    """
    Define the reaction summary.

    See Also
    --------
    Summary : Parent that defines further attributes.
    MetaboliteSummary
    ModelSummary

    """

    def __init__(
        self,
        *,
        reaction: "Reaction",
        model: "Model",
        solution: Optional["Solution"] = None,
        fva: Optional[Union[float, pd.DataFrame]] = None,
        **kwargs,
    ) -> None:
        """
        Initialize a reaction summary.

        Parameters
        ----------
        reaction: cobra.Reaction
            The reaction object whose summary we intend to get.
        model : cobra.Model
            The metabolic model in which to generate a reaction summary.
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
        ModelSummary

        """
        super(ReactionSummary, self).__init__(**kwargs)
        self._reaction = reaction.copy()
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
                model,
                reaction_list=[self._reaction],
                fraction_of_optimum=fva,
            )

        # Create the basic flux table.
        self._flux = pd.DataFrame(
            data={"flux": [solution[self._reaction.id]]},
            index=[self._reaction.id],
        )
        if fva is not None:
            self._flux = self._flux.join(fva)

    def _string_flux(
        self,
        threshold: float,
        float_format: str,
    ) -> str:
        """
        Transform a flux data frame to a string.

        Parameters
        ----------
        threshold : float
            Hide fluxes below the threshold from being displayed.
        float_format : str
            Format string for floats.

        Returns
        -------
        str
            A string representation of the flux (with ranges).

        """

        if "minimum" in self._flux.columns and "maximum" in self._flux.columns:
            frame = self._flux.loc[
                (self._flux["flux"].abs() >= threshold)
                | (self._flux["minimum"].abs() >= threshold)
                | (self._flux["maximum"].abs() >= threshold),
                :,
            ].copy()
            return (
                f"{frame.at[self._reaction.id, 'flux']:{float_format}} "
                f"[{frame.at[self._reaction.id, 'minimum']:{float_format}}; "
                f"{frame.at[self._reaction.id, 'maximum']:{float_format}}]"
            )
        else:
            frame = self._flux.loc[self._flux["flux"].abs() >= threshold, :].copy()
            return f"{frame.at[self._reaction.id, 'flux']:{float_format}}"

    def to_string(
        self,
        names: bool = False,
        threshold: Optional[float] = None,
        float_format: str = ".4G",
        column_width: int = 79,
    ) -> str:
        """
        Return a pretty string representation of the reaction summary.

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
            header = shorten(self._reaction.name, width=column_width, placeholder="...")
        else:
            header = shorten(self._reaction.id, width=column_width, placeholder="...")

        flux = self._string_flux(threshold, float_format)

        return dedent(
            f"""
            {header}
            {'=' * len(header)}
            {self._reaction.build_reaction_string(use_metabolite_names=names)}
            Bounds: {self._reaction.lower_bound}, {self._reaction.upper_bound}
            Flux: {flux}
            """
        )

    def to_html(
        self,
        names: bool = False,
        threshold: Optional[float] = None,
        float_format: str = ".4G",
    ) -> str:
        """
        Return a rich HTML representation of the reaction summary.

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
            header = self._reaction.name
        else:
            header = self._reaction.id

        flux = self._string_flux(threshold, float_format)

        return (
            f"<h3>{header}</h3>"
            f"<p>{self._reaction.build_reaction_string(use_metabolite_names=names)}</p>"
            f"<p>Bounds: {self._reaction.lower_bound}, {self._reaction.upper_bound}</p>"
            f"<p>Flux: {flux}</p>"
        )
