"""Provide a reaction summary class."""


from textwrap import dedent
from typing import TYPE_CHECKING, Optional

import pandas as pd

from cobra.flux_analysis import flux_variability_analysis
from cobra.summary import Summary


if TYPE_CHECKING:
    from cobra import Model, Reaction


class ReactionSummary(Summary):
    """
    Define the reaction summary.

    Attributes
    ----------
    flux: pandas.DataFrame
        The reaction flux (with minimum and maximum if an FVA was requested).

    See Also
    --------
    Summary : Parent that defines further attributes.
    MetaboliteSummary
    ModelSummary

    """

    def __init__(self, *, reaction: "Reaction", model: "Model", **kwargs) -> None:
        """
        Initialize a metabolite summary.

        Parameters
        ----------
        reaction: cobra.Reaction
            The reaction object whose summary we intend to get.
        model : cobra.Model
            The metabolic model in which to generate a reaction summary.

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
        super(ReactionSummary, self).__init__(model=model, **kwargs)
        self._reaction = reaction.copy()
        self.flux: Optional[pd.DataFrame] = None
        self._generate()

    def _generate(self) -> None:
        """
        Returns
        -------
        rxn_summary: pandas.DataFrame
            The DataFrame of reaction summary data.

        """
        super()._generate()

        if isinstance(self._fva, float):
            self._fva = flux_variability_analysis(
                self._model,
                reaction_list=[self._reaction],
                fraction_of_optimum=self._fva,
            )

        self.flux = pd.DataFrame(
            data={"flux": [self._solution[self._reaction.id]]},
            index=[self._reaction.id],
        )
        if self._fva is not None:
            self.flux = self.flux.join(self._fva)

    def to_string(self, names: bool = False, float_format: str = ".3G") -> str:
        if self._fva is None:
            flux = f"Flux: {self.flux.at[self._reaction.id, 'flux']:{float_format}}"
        else:
            flux = (
                f"Flux: {self.flux.at[self._reaction.id, 'flux']:{float_format}} "
                f"[{self.flux.at[self._reaction.id, 'minimum']:{float_format}}; "
                f"{self.flux.at[self._reaction.id, 'maximum']:{float_format}}]"
            )

        return dedent(
            f"""
            {flux}
            {self._reaction.build_reaction_string(use_metabolite_names=names)}
            """
        )

    def to_html(self, names: bool = False) -> str:
        return f"""
            {self.flux.to_html()}
            {self._reaction._repr_html_()}
            """
