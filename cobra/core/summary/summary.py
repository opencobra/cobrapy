# -*- coding: utf-8 -*-

"""Define the Summary class."""

from __future__ import absolute_import, division

from warnings import warn

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

    Methods
    -------
    to_frame
        Return a data frame representation of the summary.

    """

    def __init__(
            self,
            model,
            solution=None,
            threshold=None,
            fva=None,
            names=False,
            float_format="{:.3G}".format,
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
            Format string for floats (default ``'{:3G}'.format``).

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

    def _generate(self):
        """Generate the summary for the required cobra object.

        This is an abstract method and thus the subclass needs to
        implement it.

        """
        raise NotImplementedError("This method needs to be implemented by the "
                                  "subclass.")

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

        abs_flux = flux_dataframe['flux'].abs()
        flux_threshold = self.threshold * abs_flux.max()

        # Drop unused boundary fluxes i.e., fluxes below threshold
        if self.fva is None:
            flux_dataframe = \
                flux_dataframe.loc[abs_flux >= flux_threshold, :].copy()
        else:
            flux_dataframe = (
                flux_dataframe
                .loc[(abs_flux >= flux_threshold) |
                     (flux_dataframe['fmin'].abs() >= flux_threshold) |
                     (flux_dataframe['fmax'].abs() >= flux_threshold), :]
                .copy()
                )

        # Make all fluxes positive while maintaining proper direction
        if self.fva is None:
            # add information regarding direction
            flux_dataframe['is_input'] = (flux_dataframe['flux'] >= 0)
            # make the fluxes absolute
            flux_dataframe['flux'] = flux_dataframe['flux'].abs()
            # sort the values
            flux_dataframe.sort_values(by=['is_input', 'flux', 'id'],
                                       ascending=[False, False, True],
                                       inplace=True)
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

            flux_dataframe['is_input'] = (sign == 1)

            flux_dataframe.loc[:, ['flux', 'fmin', 'fmax']] = (
                flux_dataframe.loc[:, ['flux', 'fmin', 'fmax']]
                .multiply(sign, axis=0)
                .astype('float')
            )

            flux_dataframe.sort_values(by=['is_input', 'flux', 'fmax', 'fmin',
                                           'id'],
                                       ascending=[False, False, False, False,
                                                  True],
                                       inplace=True)

        return flux_dataframe

    def to_frame(self):
        """Generate a pandas DataFrame.

        This is an abstract method and thus the subclass needs to
        implement it.

        """
        raise NotImplementedError("This method needs to be implemented by the "
                                  "subclass.")

    def _to_table(self):
        """Generate a pretty-print table.

        This is an abstract method and thus the subclass needs to
        implement it.

        """
        raise NotImplementedError("This method needs to be implemented by the "
                                  "subclass.")

    def __str__(self):
        if self.float_format is not None:
            warn("Setting float_format to anything other than None "
                 "will cause nan to be present in the output.")
        return self._to_table()

    def _repr_html_(self):
        return self.to_frame()._repr_html_()
