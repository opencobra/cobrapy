"""Provide the abstract base summary class."""


import logging
from abc import ABC, abstractmethod
from typing import TYPE_CHECKING, Optional, Union

from cobra.flux_analysis import pfba
from cobra.flux_analysis.helpers import normalize_cutoff


if TYPE_CHECKING:
    from pandas import DataFrame

    from cobra import Model, Solution


logger = logging.getLogger(__name__)


class Summary(ABC):
    """
    Define the abstract base summary.

    Attributes
    ----------
    model : cobra.Model
        The metabolic model in which to generate a summary description.
    solution : cobra.Solution, optional
        A solution that matches the given model. If missing, it will be generated using
        parsimonious FBA.
    threshold : float, optional
        Threshold below which fluxes are not reported.
    fva : pandas.DataFrame or float, optional
        The result of a flux variability analysis (FVA) involving reactions of
        interest for the summary. A number will be used as fraction of optimum flux
        for calculating a new FVA.
    names : bool, optional
        Whether or not to use object names rather than identifiers.
    float_format : callable
        Format string for displaying floats.
    data_frame : pandas.DataFrame
        The table containing the summary.

    """

    def __init__(
        self,
        *,
        model: "Model",
        solution: Optional["Solution"] = None,
        threshold: Optional[float] = None,
        fva: Optional[Union[float, "DataFrame"]] = None,
        **kwargs,
    ) -> None:
        """
        Initialize a summary.

        Parameters
        ----------
        model : cobra.Model
            The metabolic model in which to generate a summary description.
        solution : cobra.Solution, optional
            A previous model solution to use for generating the summary. If
            None, the summary method will resolve the model.  Note that the
            solution object must match the model, i.e., changes to the model
            such as changed bounds, added or removed reactions are not taken
            into account by this method (default None).
        threshold : float, optional
            Threshold below which fluxes are not reported. May not be smaller
            than the model tolerance (default None).
        fva : pandas.DataFrame or float, optional
            Whether or not to include flux variability analysis in the output.
            If given, fva should either be a previous FVA solution matching the
            model or a float between 0 and 1 representing the fraction of the
            optimum objective to be searched (default None).
        names : bool, optional
            Emit reaction and metabolite names rather than identifiers (default
            False).
        float_format : callable, optional
            Format string for floats (default ``'{:3G}'.format``). Please see
            https://pandas.pydata.org/pandas-docs/stable/user_guide/style.html#Finer-Control:-Display-Values
            for more information.

        Other Parameters
        ----------------
        kwargs :
            Further keyword arguments are passed on to the parent class.

        """
        super().__init__(**kwargs)
        self._model = model
        self._solution = solution
        self._threshold = normalize_cutoff(self._model, threshold)
        self._fva = fva
        self._flux_frame = None

    def _generate(self) -> None:
        """Generate the summary for the required cobra object.

        This is an abstract method and thus the subclass needs to
        implement it.

        """
        if self._solution is None:
            logger.info("Generating new parsimonious flux distribution.")
            self._solution = pfba(self._model)

    @abstractmethod
    def to_frame(self, names: bool = False):
        """Return a pandas data frame representation of the summary."""
        raise NotImplementedError(
            "This method needs to be implemented by the subclass."
        )

    @abstractmethod
    def to_html(self, names: bool = False):
        """Return a rich HTML representation of the summary."""
        raise NotImplementedError(
            "This method needs to be implemented by the subclass."
        )

    @abstractmethod
    def to_string(self, names: bool = False):
        """Return a string representation of the summary."""
        raise NotImplementedError(
            "This method needs to be implemented by the subclass."
        )

    def __str__(self) -> str:
        return self.to_string()

    def _repr_html_(self) -> str:
        return self.to_html()

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
                """Decide whether or not to reverse a flux to make it
                positive."""

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

    # def __str__(self):
    #     return self._display().to_string(
    #         header=True, index=True, na_rep="", formatters=self._styles()
    #     )

    # def _repr_html_(self):
    #     self._display().style.format(self._styles())
