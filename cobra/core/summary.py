# -*- coding: utf-8 -*-

"""Define the Summary class."""

from __future__ import absolute_import, division

import logging
from operator import attrgetter

import numpy as np
import pandas as pd
from six import iteritems, print_
from six.moves import zip_longest
from tabulate import tabulate

from cobra.core import get_solution
from cobra.flux_analysis.variability import flux_variability_analysis
from cobra.util.solver import linear_reaction_coefficients
from cobra.util.util import format_long_string


LOGGER = logging.getLogger(__name__)


class Summary(object):
    """Class definition for a Summary object.

    Parameters
    ----------
    solution : cobra.Solution or None
        A previously solved model solution to use for generating the
        summary. If None, the summary method will resolve the model.
        Note that the solution object must match the model, i.e., changes
        to the model such as changed bounds, added or removed reactions are
        not taken into account by this method.
    threshold : float
        Threshold below which fluxes are not reported.
    fva : pandas.DataFrame, float or None
        Whether or not to include flux variability analysis in the output.
        If given, fva should either be a previous FVA solution matching
        the model or a float between 0 and 1 representing the
        fraction of the optimum objective to be searched.
    names : bool
        Emit reaction and metabolite names rather than identifiers.
    float_format : one-parameter function
        Format string for floats.

    """

    def __init__(self, solution, threshold, fva, names, float_format):
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
            flux_dataframe['fva_fmt'] = \
                flux_dataframe.apply(lambda x: ("[{0.fmin:" + self.floatfmt + "},\
                {0.fmax:" + self.floatfmt + "}]").format(x), axis=1)

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


class MetaboliteSummary(Summary):
    """Class definition for a MetaboliteSummary object.

    Parameters
    ----------
    met: cobra.Metabolite
        The Metabolite object whose summary we intend to get.
    solution : cobra.Solution or None
        A previously solved model solution to use for generating the
        summary. If None, the summary method will resolve the model.
        Note that the solution object must match the model, i.e., changes
        to the model such as changed bounds, added or removed reactions are
        not taken into account by this method.
    threshold : float
        Threshold below which fluxes are not reported.
    fva : pandas.DataFrame, float or None
        Whether or not to include flux variability analysis in the output.
        If given, fva should either be a previous FVA solution matching
        the model or a float between 0 and 1 representing the
        fraction of the optimum objective to be searched.
    names : bool
        Emit reaction and metabolite names rather than identifiers.
    float_format : one-parameter function
        Format string for floats.

    """

    def __init__(self, met, solution, threshold, fva, names, float_format):
        super(MetaboliteSummary, self).__init__(solution, threshold, fva,
                                                names, float_format)
        self.met = met

    def _generate(self):
        """
        Returns
        -------
        flux_summary: pandas.DataFrame
            The DataFrame of flux summary data.

        """
        if self.names:
            emit = attrgetter('name')
        else:
            emit = attrgetter('id')

        if self.solution is None:
            self.met.model.slim_optimize(error_value=None)
            self.solution = get_solution(self.met.model,
                                         reactions=self.met.reactions)

        rxns = sorted(self.met.reactions, key=attrgetter("id"))
        rxn_id = list()
        rxn_name = list()
        flux = list()
        reaction = list()

        for rxn in rxns:
            rxn_id.append(rxn.id)
            rxn_name.append(format_long_string(emit(rxn), 10))
            flux.append(self.solution[rxn.id] * rxn.metabolites[self.met])
            txt = rxn.build_reaction_string(use_metabolite_names=self.names)
            reaction.append(format_long_string(txt, 40
                                               if self.fva is not None
                                               else 50))

        flux_summary = pd.DataFrame({
            "id": rxn_name,
            "flux": flux,
            "reaction": reaction
        }, index=rxn_id)

        if self.fva is not None:
            if hasattr(self.fva, 'columns'):
                fva_results = self.fva
            else:
                fva_results = flux_variability_analysis(
                    self.met.model, list(self.met.reactions),
                    fraction_of_optimum=self.fva)

            flux_summary["maximum"] = np.zeros(len(rxn_id), dtype=float)
            flux_summary["minimum"] = np.zeros(len(rxn_id), dtype=float)

            for rxn in rxns:
                fmax = rxn.metabolites[self.met] * \
                    fva_results.at[rxn.id, "maximum"]

                fmin = rxn.metabolites[self.met] * \
                    fva_results.at[rxn.id, "minimum"]

                if abs(fmin) <= abs(fmax):
                    flux_summary.at[rxn.id, "fmax"] = fmax
                    flux_summary.at[rxn.id, "fmin"] = fmin
                else:
                    # Reverse fluxes.
                    flux_summary.at[rxn.id, "fmax"] = fmin
                    flux_summary.at[rxn.id, "fmin"] = fmax

        assert flux_summary["flux"].sum() < 1E-6, "Error in flux balance"

        flux_summary = self._process_flux_dataframe(flux_summary)

        flux_summary['percent'] = 0
        total_flux = flux_summary.loc[flux_summary.is_input, "flux"].sum()

        flux_summary.loc[flux_summary.is_input, 'percent'] = \
            flux_summary.loc[flux_summary.is_input, 'flux'] / total_flux

        flux_summary.loc[~flux_summary.is_input, 'percent'] = \
            flux_summary.loc[~flux_summary.is_input, 'flux'] / total_flux

        flux_summary['percent'] = flux_summary.percent.\
            apply(lambda x: '{:.0%}'.format(x))

        return flux_summary

    def to_frame(self):
        """
        Returns
        -------
        A pandas.DataFrame of the summary.

        """
        flux_df = self._generate()

        if self.fva is not None:
            column_names = ['percent', 'flux', 'fmin',
                            'fmax', 'id', 'reaction']

            # Generate DataFrame of production reactions
            flux_prod = flux_df[flux_df.is_input.values]\
                .loc[:, column_names]\
                .reset_index(drop=True)
            flux_prod.columns = [name.upper() for name in column_names]
            flux_prod['RXN_STAT'] = 'PRODUCING'

            # Generate DataFrame of consumption reactions
            flux_cons = flux_df[~flux_df.is_input.values]\
                .loc[:, column_names]\
                .reset_index(drop=True)
            flux_cons.columns = [name.upper() for name in column_names]
            flux_cons['RXN_STAT'] = 'CONSUMING'

            concat_df = pd.concat([flux_prod, flux_cons])

            del flux_prod, flux_cons

            concat_df.rename(columns={'FMIN': 'FLUX_MIN', 'FMAX': 'FLUX_MAX'},
                             inplace=True)
            concat_df.set_index(['RXN_STAT', 'ID'], inplace=True)

        else:
            column_names = ['percent', 'flux', 'id', 'reaction']

            # Generate DataFrame of production reactions
            flux_prod = flux_df[flux_df.is_input.values]\
                .loc[:, column_names]\
                .reset_index(drop=True)
            flux_prod.columns = [name.upper() for name in column_names]
            flux_prod['RXN_STAT'] = 'PRODUCING'

            # Generate DataFrame of consumption reactions
            flux_cons = flux_df[~flux_df.is_input.values]\
                .loc[:, column_names]\
                .reset_index(drop=True)
            flux_cons.columns = [name.upper() for name in column_names]
            flux_cons['RXN_STAT'] = 'CONSUMING'

            concat_df = pd.concat([flux_prod, flux_cons])

            del flux_prod, flux_cons

            concat_df.set_index(['RXN_STAT', 'ID'], inplace=True)

        return concat_df

    def _to_table(self):
        """
        Returns
        -------
        A string of the summary table.

        """
        return self.to_frame().to_string(header=True, index=True, na_rep='',
                                         float_format=self.float_format,
                                         sparsify=True, justify='center')


class ModelSummary(Summary):
    """Class definition for a ModelSummary object.

    Parameters
    ----------
    model: cobra.Model
        The Model object whose summary we intend to get.
    solution : cobra.Solution or None
        A previously solved model solution to use for generating the
        summary. If None, the summary method will resolve the model.
        Note that the solution object must match the model, i.e., changes
        to the model such as changed bounds, added or removed reactions are
        not taken into account by this method.
    threshold : float
        Threshold below which fluxes are not reported.
    fva : pandas.DataFrame, float or None
        Whether or not to include flux variability analysis in the output.
        If given, fva should either be a previous FVA solution matching
        the model or a float between 0 and 1 representing the
        fraction of the optimum objective to be searched.
    names : bool
        Emit reaction and metabolite names rather than identifiers.
    float_format : one-parameter function
        Format string for floats.

    """

    def __init__(self, model, solution, threshold, fva, names, float_format):
        super(ModelSummary, self).__init__(solution, threshold, fva, names,
                                           float_format)
        self.model = model

    def _generate(self):
        """
        Returns
        -------
        metabolite_fluxes: pandas.DataFrame
            The DataFrame of metabolite fluxes.
        obj_fluxes: pandas.DataFrame
            The DataFrame of objective fluxes.

        """
        if self.names:
            emit = attrgetter('name')
        else:
            emit = attrgetter('id')

        objective_reactions = linear_reaction_coefficients(self.model)
        boundary_reactions = self.model.exchanges
        summary_rxns = \
            set(objective_reactions.keys()).union(boundary_reactions)

        if self.solution is None:
            self.model.slim_optimize(error_value=None)
            self.solution = get_solution(self.model, reactions=summary_rxns)

        # Create a pandas.DataFrame of objective fluxes
        obj_fluxes = pd.DataFrame({key: self.solution[key.id] * value
                                   for key, value in
                                   iteritems(objective_reactions)},
                                  index=['flux']).T
        obj_fluxes['id'] = obj_fluxes\
            .apply(lambda x: format_long_string(x.name.id, 15), axis=1)

        # Build a dictionary of metabolite production from the
        # boundary reactions
        metabolites = {met for rxn in boundary_reactions
                       for met in rxn.metabolites}
        index = sorted(metabolites, key=attrgetter('id'))

        # Create a pandas.DataFrame of metabolite fluxes
        metabolite_fluxes = pd.DataFrame({
            'id': [format_long_string(emit(met), 15) for met in index],
            'flux': np.zeros(len(index), dtype=float)
        }, index=[met.id for met in index])

        # Set the proper flux values for metabolites
        for rxn in boundary_reactions:
            for met, stoich in iteritems(rxn.metabolites):
                metabolite_fluxes.at[met.id, 'flux'] += \
                    stoich * self.solution[rxn.id]

        # Calculate FVA results if requested
        if self.fva is not None:
            if len(index) != len(boundary_reactions):
                LOGGER.warning(
                    "There exists more than one boundary reaction per "
                    "metabolite. Please be careful when evaluating flux "
                    "ranges.")
            metabolite_fluxes['fmin'] = np.zeros(len(index), dtype=float)
            metabolite_fluxes['fmax'] = np.zeros(len(index), dtype=float)

            if hasattr(self.fva, 'columns'):
                fva_results = self.fva
            else:
                fva_results = flux_variability_analysis(
                    self.model, reaction_list=boundary_reactions,
                    fraction_of_optimum=self.fva)

            for rxn in boundary_reactions:
                for met, stoich in iteritems(rxn.metabolites):
                    fmin = stoich * fva_results.at[rxn.id, 'minimum']
                    fmax = stoich * fva_results.at[rxn.id, 'maximum']
                    # Correct 'max' and 'min' for negative values
                    if abs(fmin) <= abs(fmax):
                        metabolite_fluxes.at[met.id, 'fmin'] += fmin
                        metabolite_fluxes.at[met.id, 'fmax'] += fmax
                    else:
                        metabolite_fluxes.at[met.id, 'fmin'] += fmax
                        metabolite_fluxes.at[met.id, 'fmax'] += fmin

        # Generate a pandas.DataFrame of boundary fluxes
        metabolite_fluxes = self._process_flux_dataframe(metabolite_fluxes)

        return metabolite_fluxes, obj_fluxes

    def to_frame(self):
        """
        Returns
        -------
        A pandas.DataFrame of the summary.

        """
        met_df, obj_df = self._generate()

        if self.fva is not None:
            column_names = [['IN_FLUXES', 'OUT_FLUXES', 'OBJECTIVES'],
                            ['ID', 'FLUX', 'FLUX_MIN', 'FLUX_MAX']]

            # Generate new DataFrames for easy concatenation
            met_in_df = met_df[met_df['is_input']]\
                .loc[:, ['id', 'flux', 'fmin', 'fmax']].reset_index(drop=True)

            met_out_df = met_df[~met_df['is_input']]\
                .loc[:, ['id', 'flux', 'fmin', 'fmax']].reset_index(drop=True)

            obj_df = obj_df.loc[:, ['id', 'flux']].reset_index(drop=True)

            # concatenate DataFrames
            concat_df = pd.concat([met_in_df, met_out_df, obj_df], axis=1)\
                          .values

            del met_in_df, met_out_df, obj_df

            # Generate column names
            columns = pd.MultiIndex.from_product(column_names)

            # Remove 'FLUX_MIN' and 'FLUX_MAX' for 'OBJECTIVES'
            columns.set_codes([[0, 0, 0, 0, 2, 2, 2, 2, 1, 1],
                               [3, 0, 2, 1, 3, 0, 2, 1, 3, 0]],
                              inplace=True,
                              verify_integrity=False)

        else:
            column_names = [['IN_FLUXES', 'OUT_FLUXES', 'OBJECTIVES'],
                            ['ID', 'FLUX']]

            # Generate new DataFrames for easy concatenation
            met_in_df = met_df[met_df['is_input']]\
                .loc[:, ['id', 'flux']].reset_index(drop=True)

            met_out_df = met_df[~met_df['is_input']]\
                .loc[:, ['id', 'flux']].reset_index(drop=True)

            obj_df = obj_df.loc[:, ['id', 'flux']].reset_index(drop=True)

            # concatenate DataFrames
            concat_df = pd.concat([met_in_df, met_out_df, obj_df], axis=1)\
                          .values

            del met_in_df, met_out_df, obj_df

            # Generate column names
            columns = pd.MultiIndex.from_product(column_names)

        return pd.DataFrame(
            data=concat_df,
            columns=columns
        )

    def _to_table(self):
        """
        Returns
        -------
        A string of the summary table.

        """
        return self.to_frame().to_string(header=True, index=False, na_rep='',
                                         float_format=self.float_format,
                                         sparsify=False, justify='center')


class ReactionSummary(Summary):
    """Class definition for a ReactionSummary object.

    Parameters
    ----------
    rxn: cobra.Reaction
        The Reaction object whose summary we intend to get.
    names : bool, optional
        Emit gene and metabolite names rather than identifiers (default
        False).

    """

    def __init__(self, rxn, names):
        super(ReactionSummary, self).__init__(names=names, solution=None,
                                              threshold=None, fva=None,
                                              floatfmt=None)
        self.rxn = rxn

    def _generate(self):
        """
        Returns
        -------
        rxn_summary: pandas.DataFrame
            The DataFrame of reaction summary data.

        """
        if self.names:
            emit = attrgetter('name')
        else:
            emit = attrgetter('id')

        gene_temp_df = pd.DataFrame([emit(gene) for gene in self.rxn.genes])
        met_temp_df = pd.DataFrame([
            [emit(key), value, key.compartment] for key, value in
            iteritems(self.rxn.metabolites)
        ])
        data = pd.concat([gene_temp_df, met_temp_df], axis=1).fillna('').values

        columns = pd.MultiIndex.from_tuples((('REACTION', 'GENES', 'ID'),
                                             ('REACTION', 'METABOLITES', 'ID'),
                                             ('REACTION', 'METABOLITES',
                                              'STOICHIOMETRIC COEFFICIENT'),
                                             ('REACTION', 'METABOLITES',
                                              'COMPARTMENT')))

        rxn_summary = pd.DataFrame(
            data=data,
            columns=columns
        )

        del gene_temp_df, met_temp_df

        return rxn_summary

    def to_table(self):
        """
        Returns
        -------
        Nothing

        """
        rxn_df = self._generate()

        gene_table = tabulate(
            rxn_df['REACTION', 'GENES'].replace('', np.nan).dropna().values,
            headers=['ID']
        )

        reactants_table = tabulate(
            rxn_df[rxn_df['REACTION', 'METABOLITES',
                          'STOICHIOMETRIC COEFFICIENT'] < 0]
            .loc[:, ('REACTION', 'METABOLITES')].values,
            headers=['ID', 'STOICHIOMETRIC COEFFICIENT', 'COMPARTMENT']
        )

        products_table = tabulate(
            rxn_df[rxn_df['REACTION', 'METABOLITES',
                          'STOICHIOMETRIC COEFFICIENT'] > 0]
            .loc[:, ('REACTION', 'METABOLITES')].values,
            headers=['ID', 'STOICHIOMETRIC COEFFICIENT', 'COMPARTMENT']
        )

        rxn_tag = '{0} {1}'.format(format_long_string(self.rxn.name, 45),
                                   format_long_string(self.rxn.id, 10))

        head = 'REACTANTS -- ' + rxn_tag

        print_('REACTION: ' +
               self.rxn.build_reaction_string(use_metabolite_names=self.names))

        print_()
        print_('GENES -- ' + rxn_tag)
        print_('-' * len(head))
        print_(gene_table)

        print_()
        print_(head)
        print_('-' * len(head))
        print_(reactants_table)

        print_()
        print_('PRODUCTS -- ' + rxn_tag)
        print_('-' * len(head))
        print_(products_table)

    def to_frame(self):
        """
        Returns
        -------
        A pandas.DataFrame of the summary.

        """
        return self._generate()
