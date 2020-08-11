"""Provide the abstract base summary class."""


import logging
from abc import ABC, abstractmethod
from typing import TYPE_CHECKING, Optional, Union


if TYPE_CHECKING:
    from pandas import DataFrame

    from cobra import Model, Solution


logger = logging.getLogger(__name__)


class Summary(ABC):
    """
    Define the abstract base summary.

    See Also
    --------
    MetaboliteSummary
    ReactionSummary
    ModelSummary

    """

    def __init__(self, **kwargs,) -> None:
        """
        Initialize a summary.

        Other Parameters
        ----------------
        kwargs :
            Further keyword arguments are passed on to the parent class.

        """
        super().__init__(**kwargs)
        self._flux = None

    @abstractmethod
    def _generate(
        self,
        model: "Model",
        solution: Optional["Solution"],
        fva: Optional[Union[float, "DataFrame"]],
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
        raise NotImplementedError(
            "This method needs to be implemented by the subclass."
        )

    def __str__(self) -> str:
        """Return a string representation of the summary."""
        return self.to_string()

    def _repr_html_(self) -> str:
        """Return a rich HTML representation of the summary."""
        return self.to_html()

    @abstractmethod
    def to_string(
        self,
        names: bool = False,
        threshold: float = 1e-6,
        float_format: str = ".4G",
        column_width: int = 79,
    ) -> str:
        """
        Return a pretty string representation of the summary.

        Parameters
        ----------
        names : bool, optional
            Whether or not elements should be displayed by their common names
            (default False).
        threshold : float, optional
            Hide fluxes below the threshold from being displayed (default 1e-6).
        float_format : str, optional
            Format string for floats (default '.4G').
        column_width : int, optional
            The maximum column width for each row (default 79).

        Returns
        -------
        str
            The summary formatted as a pretty string.

        """
        raise NotImplementedError(
            "This method needs to be implemented by the subclass."
        )

    @abstractmethod
    def to_html(
        self, names: bool = False, threshold: float = 1e-6, float_format: str = ".4G"
    ) -> str:
        """
        Return a rich HTML representation of the metabolite summary.

        Parameters
        ----------
        names : bool, optional
            Whether or not elements should be displayed by their common names
            (default False).
        threshold : float, optional
            Hide fluxes below the threshold from being displayed (default 1e-6).
        float_format : str, optional
            Format string for floats (default '.4G').

        Returns
        -------
        str
            The summary formatted as HTML.

        """
        raise NotImplementedError(
            "This method needs to be implemented by the subclass."
        )

    def to_frame(self) -> "DataFrame":
        """Return the a data frame representation of the summary."""
        return self._flux.copy()

    def _process_flux_dataframe(self, flux_dataframe):
        """Process a flux DataFrame for convenient downstream analysis.

        This method removes flux entries which are below the threshold and
        also adds information regarding the direction of the fluxes. It is used
        in both ModelSummary and MetaboliteSummary.

        Parameters
        ----------
        flux_dataframe: pandas.DataFrame
            The pandas.DataFrame to process.

        Returns
        -------
        A processed pandas.DataFrame.

        """

        abs_flux = flux_dataframe["flux"].abs()
        flux_threshold = self.threshold * abs_flux.max()

        # Drop unused boundary fluxes i.e., fluxes below threshold
        if self.fva is None:
            flux_dataframe = flux_dataframe.loc[abs_flux >= flux_threshold, :].copy()
        else:
            flux_dataframe = flux_dataframe.loc[
                (abs_flux >= flux_threshold)
                | (flux_dataframe["fmin"].abs() >= flux_threshold)
                | (flux_dataframe["fmax"].abs() >= flux_threshold),
                :,
            ].copy()

        # Make all fluxes positive while maintaining proper direction
        if self.fva is None:
            # add information regarding direction
            flux_dataframe["is_input"] = flux_dataframe["flux"] >= 0
            # make the fluxes absolute
            flux_dataframe["flux"] = flux_dataframe["flux"].abs()
            # sort the values
            flux_dataframe.sort_values(
                by=["is_input", "flux", "id"],
                ascending=[False, False, True],
                inplace=True,
            )
        else:

            def get_direction(row):
                """Decide whether or not to reverse a flux to make it positive."""

                if row.flux < 0:
                    sign = -1
                elif row.flux > 0:
                    sign = 1
                elif (row.fmax > 0) & (row.fmin <= 0):
                    sign = 1
                elif (row.fmax < 0) & (row.fmin >= 0):
                    sign = -1
                elif ((row.fmax + row.fmin) / 2) < 0:
                    sign = -1
                else:
                    sign = 1

                return sign

            # get a sign DataFrame to use as a pseudo-mask
            sign = flux_dataframe.apply(get_direction, axis=1)

            flux_dataframe["is_input"] = sign == 1

            flux_dataframe.loc[:, ["flux", "fmin", "fmax"]] = (
                flux_dataframe.loc[:, ["flux", "fmin", "fmax"]]
                .multiply(sign, axis=0)
                .astype("float")
            )

            flux_dataframe.sort_values(
                by=["is_input", "flux", "fmax", "fmin", "id"],
                ascending=[False, False, False, False, True],
                inplace=True,
            )

        return flux_dataframe
