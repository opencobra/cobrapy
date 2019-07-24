# -*- coding: utf-8 -*-

"""Define the Summary class."""

from __future__ import absolute_import, division

from cobra.flux_analysis import pfba
from cobra.flux_analysis.helpers import normalize_cutoff


class Summary(object):
    """
    Define the abstract base summary.

    Attributes
    ----------
    model : cobra.Model
        The metabolic model in which to generate a summary description.
    solution : cobra.Solution
        A solution that matches the given model.
    threshold : float, optional
        Threshold below which fluxes are not reported.
    fva : pandas.DataFrame, optional
        The result of a flux variability analysis (FVA) involving reactions of
        interest if an FVA was requested.
    names : bool
        Whether or not to use object names rather than identifiers.
    float_format : callable
        Format string for displaying floats.
    data_frame : pandas.DataFrame
        The table containing the summary.

    """

    def __init__(
            self,
            model,
            solution=None,
            threshold=None,
            fva=None,
            names=False,
            float_format='{:.3G}'.format,
            **kwargs
    ):
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
        super(Summary, self).__init__(**kwargs)
        self.model = model
        self.solution = solution
        self.threshold = normalize_cutoff(self.model, threshold)
        self.fva = fva
        self.names = names
        self.float_format = float_format
        self.data_frame = None
        self._generate()

    def __repr__(self):
        return "<%s in %s>".format(type(self).__name__, str(self.model))

    def __str__(self):
        return self._display().to_string(
            header=True,
            index=True,
            na_rep='',
            formatters=self._styles()
        )

    def _repr_html_(self):
        self._display().style.format(self._styles())

    def _display(self):
        raise NotImplementedError("Abstract base method.")

    def _styles(self):
        """Set the display options on the data frame."""
        styles = {'flux': self.float_format, 'percent': '{:.2%}'}
        if self.fva is not None:
            styles['minimum'] = self.float_format
            styles['maximum'] = self.float_format
        return styles

    def _generate(self):
        """Generate the summary."""
        if self.solution is None:
            self.solution = pfba(self.model)
