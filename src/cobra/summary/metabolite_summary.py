# -*- coding: utf-8 -*-

"""Define the MetaboliteSummary class."""

from __future__ import absolute_import, division

from operator import attrgetter

from pandas import DataFrame

from cobra.flux_analysis.variability import flux_variability_analysis
from cobra.summary import Summary


class MetaboliteSummary(Summary):
    """
    Define the metabolite summary.

    Attributes
    ----------
    metabolite: cobra.Metabolite
        The metabolite to summarize.

    See Also
    --------
    Summary : Parent that defines further attributes.
    ReactionSummary
    ModelSummary

    """

    def __init__(self, metabolite, model, **kwargs):
        """
        Initialize a metabolite summary.

        Parameters
        ----------
        metabolite: cobra.Metabolite
            The metabolite object whose summary we intend to get.
        model : cobra.Model
            The metabolic model for which to generate a metabolite summary.

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
        super(MetaboliteSummary, self).__init__(model=model, **kwargs)
        self.metabolite = metabolite
        self._generate()

    def _display(self):
        output = self.data_frame.copy()
        if self.names:
            output["name"] = [
                self.model.reactions.get_by_id(r).name for r in output.index
            ]
        output["definition"] = [
            self.model.reactions.get_by_id(r).build_reaction_string(self.names)
            for r in output.index
        ]
        return output

    def _generate(self):
        """
        Returns
        -------
        flux_summary: pandas.DataFrame
            The DataFrame of flux summary data.

        """
        rxns = sorted(self.metabolite.reactions, key=attrgetter("id"))
        data = [
            (
                r.id,
                self.solution[r.id] * r.metabolites[self.metabolite],
                r.metabolites[self.metabolite],
            )
            for r in rxns
        ]
        flux_summary = DataFrame(data=data, columns=["reaction", "flux", "factor"])

        if self.fva is not None:
            if not isinstance(self.fva, DataFrame):
                self.fva = flux_variability_analysis(
                    self.model, reaction_list=rxns, fraction_of_optimum=self.fva
                )

            flux_summary = flux_summary.set_index("reaction").join(self.fva)
            flux_summary.loc[
                (
                    flux_summary[["flux", "minimum", "maximum"]].abs()
                    < self.model.tolerance
                ),
                ["flux", "minimum", "maximum"],
            ] = 0
            flux_summary[["minimum", "maximum"]] = flux_summary[
                ["minimum", "maximum"]
            ].mul(flux_summary["factor"], axis=0)
            # Negative factors inverse the inequality relationship.
            negative = flux_summary["factor"] < 0
            tmp = flux_summary.loc[negative, "maximum"]
            flux_summary.loc[negative, "maximum"] = flux_summary.loc[
                negative, "minimum"
            ]
            flux_summary.loc[negative, "minimum"] = tmp
            flux_summary = flux_summary[
                (flux_summary["flux"].abs() >= self.threshold)
                | (flux_summary["minimum"].abs() >= self.threshold)
                | (flux_summary["maximum"].abs() >= self.threshold)
            ].copy()
        else:
            flux_summary = flux_summary[
                flux_summary["flux"].abs() >= self.threshold
            ].copy()

        is_input = flux_summary["flux"] > 0
        total_flux = flux_summary.loc[is_input, "flux"].sum()
        flux_summary["percent"] = flux_summary["flux"].abs() / total_flux
        self.data_frame = flux_summary[["percent", "flux", "minimum", "maximum"]]
