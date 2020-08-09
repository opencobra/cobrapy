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
    reaction: cobra.Reaction
        The reaction to summarize.

    See Also
    --------
    Summary : Parent that defines further attributes.
    MetaboliteSummary
    ModelSummary

    """

    def __init__(self, *, reaction: "Reaction", model: "Model", **kwargs):
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
        self._flux_frame: Optional[pd.DataFrame] = None
        self._generate()

    def _generate(self):
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

        self._flux_frame = pd.DataFrame(
            data={"flux": [self._solution[self._reaction.id]]},
            index=[self._reaction.id],
        )
        if self._fva is not None:
            self._flux_frame = self._flux_frame.join(self._fva)

    def __str__(self):
        return dedent(
            f"""
            {self._reaction.build_reaction_string(use_metabolite_names=self._names)}
            """
        )

    def _repr_html_(self):
        return f"""
            {self._flux_frame.to_html()}
            {self._reaction._repr_html_()}
            """

    def to_frame(self):
        """
        Returns
        -------
        A pandas.DataFrame of the summary.

        """
        return self._flux_frame

    def _to_table(self):
        """
        Returns
        -------
        A string of the summary table.

        """
        return self.to_frame().to_string(
            header=True, index=False, na_rep="", sparsify=False, justify="center"
        )
