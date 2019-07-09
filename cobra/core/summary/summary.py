# -*- coding: utf-8 -*-

"""Define the Summary class."""

from __future__ import absolute_import, division


class Summary(object):
    """Define the abstract base summary class.

    Parameters
    ----------
    solution : cobra.Solution
        A previously solved model solution to use for generating the
        summary. If None, the summary method will resolve the model.
        Note that the solution object must match the model, i.e., changes
        to the model such as changed bounds, added or removed reactions are
        not taken into account by this method.
    threshold : float
        Threshold below which fluxes are not reported.
    fva : pandas.DataFrame, float
        Whether or not to include flux variability analysis in the output.
        If given, fva should either be a previous FVA solution matching
        the model or a float between 0 and 1 representing the
        fraction of the optimum objective to be searched.
    names : bool
        Emit reaction and metabolite names rather than identifiers.
    float_format : one-parameter function
        Format string for floats.

    """

    def __init__(self, solution, threshold, fva, names, float_format,
                 **kwargs):
        super(Summary, self).__init__(**kwargs)
        self.solution = solution
        self.threshold = threshold
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
        """Some common methods for processing a database of flux information
        into print-ready formats. Used in both ModelSummary and
        MetaboliteSummary."""

        abs_flux = flux_dataframe['flux'].abs()
        flux_threshold = self.threshold * abs_flux.max()

        # Drop unused boundary fluxes
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

            # Why set to zero? If included show true value?
            # flux_dataframe.loc[
            #     flux_dataframe['flux'].abs() < flux_threshold, 'flux'] = 0

        # Make all fluxes positive
        if self.fva is None:
            flux_dataframe['is_input'] = (flux_dataframe['flux'] >= 0)
            flux_dataframe['flux'] = flux_dataframe['flux'].abs()
        else:

            def get_direction(flux, fmin, fmax):
                """Decide whether or not to reverse a flux to make it
                positive."""

                if flux < 0:
                    return -1
                elif flux > 0:
                    return 1
                elif (fmax > 0) & (fmin <= 0):
                    return 1
                elif (fmax < 0) & (fmin >= 0):
                    return -1
                elif ((fmax + fmin) / 2) < 0:
                    return -1
                else:
                    return 1

            sign = flux_dataframe\
                .apply(lambda x: get_direction(x.flux, x.fmin, x.fmax), axis=1)

            flux_dataframe['is_input'] = sign == 1

            flux_dataframe.loc[:, ['flux', 'fmin', 'fmax']] = (
                flux_dataframe.loc[:, ['flux', 'fmin', 'fmax']]
                .multiply(sign, axis=0)
                .astype('float')
                .round(6)
                )

            flux_dataframe.loc[:, ['flux', 'fmin', 'fmax']] = (
                flux_dataframe.loc[:, ['flux', 'fmin', 'fmax']]
                .applymap(lambda x: x if abs(x) > 1E-6 else 0)
                )

        if self.fva is not None:
            flux_dataframe = (
                flux_dataframe
                .sort_values(by=['flux', 'fmax', 'fmin', 'id'],
                             ascending=[False, False, False, True])
                )

        else:
            flux_dataframe = flux_dataframe\
                .sort_values(by=['flux', 'id'], ascending=[False, True])

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
        return self._to_table()

    def _repr_html_(self):
        return self.to_frame()._repr_html_()
