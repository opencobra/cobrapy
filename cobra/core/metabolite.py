# -*- coding: utf-8 -*-

from __future__ import absolute_import

import re
from warnings import warn

from six import iteritems

from cobra.core.species import Species
from cobra.core.formula import elements_and_molecular_weights

# Numbers are not required because of the |(?=[A-Z])? block. See the
# discussion in https://github.com/opencobra/cobrapy/issues/128 for
# more details.
element_re = re.compile("([A-Z][a-z]?)([0-9.]+[0-9.]?|(?=[A-Z])?)")


class Metabolite(Species):
    """Metabolite is a class for holding information regarding
    a metabolite in a cobra.Reaction object.

    Parameters
    ----------
    id : str
        the identifier to associate with the metabolite
    formula : str
        Chemical formula (e.g. H2O)
    name : str
        A human readable name.
    charge : float
       The charge number of the metabolite
    compartment: str or None
       Compartment of the metabolite.
    """

    def __init__(self, id=None, formula=None, name="",
                 charge=None, compartment=None):
        Species.__init__(self, id, name)
        self.formula = formula
        # because in a Model a metabolite may participate in multiple Reactions
        self.compartment = compartment
        self.charge = charge

        self._constraint_sense = 'E'
        self._bound = 0.

    def _set_id_with_model(self, value):
        if value in self.model.metabolites:
            raise ValueError("The model already contains a metabolite with "
                             "the id:", value)
        self.model.solver.constraints[self.id].name = value
        self._id = value
        self.model.metabolites._generate_index()

    @property
    def constraint(self):
        """Get the constraints associated with this metabolite from the solve

        Returns
        -------
        optlang.<interface>.Constraint
            the optlang constraint for this metabolite
        """
        if self.model is not None:
            return self.model.solver.constraints[self.id]

    @property
    def elements(self):
        """ Dictionary of elements as keys and their count in the metabolite
        as integer. When set, the `formula` property is update accordingly """
        tmp_formula = self.formula
        if tmp_formula is None:
            return {}
        # necessary for some old pickles which use the deprecated
        # Formula class
        tmp_formula = str(self.formula)
        # commonly occuring characters in incorrectly constructed formulas
        if "*" in tmp_formula:
            warn("invalid character '*' found in formula '%s'" % self.formula)
            tmp_formula = tmp_formula.replace("*", "")
        if "(" in tmp_formula or ")" in tmp_formula:
            warn("invalid formula (has parenthesis) in '%s'" % self.formula)
            return None
        composition = {}
        parsed = element_re.findall(tmp_formula)
        for (element, count) in parsed:
            if count == '':
                count = 1
            else:
                try:
                    count = float(count)
                    int_count = int(count)
                    if count == int_count:
                        count = int_count
                    else:
                        warn("%s is not an integer (in formula %s)" %
                             (count, self.formula))
                except ValueError:
                    warn("failed to parse %s (in formula %s)" %
                         (count, self.formula))
                    return None
            if element in composition:
                composition[element] += count
            else:
                composition[element] = count
        return composition

    @elements.setter
    def elements(self, elements_dict):
        def stringify(element, number):
            return element if number == 1 else element + str(number)

        self.formula = ''.join(stringify(e, n) for e, n in
                               sorted(iteritems(elements_dict)))

    @property
    def formula_weight(self):
        """Calculate the formula weight"""
        try:
            return sum([count * elements_and_molecular_weights[element]
                        for element, count in self.elements.items()])
        except KeyError as e:
            warn("The element %s does not appear in the peridic table" % e)

    @property
    def y(self):
        """The shadow price for the metabolite in the most recent solution

        Shadow prices are computed from the dual values of the bounds in
        the solution.

        """
        warn("use metabolite.shadow_price instead", DeprecationWarning)
        try:
            return self._model.solution.y_dict[self.id]
        except Exception as e:
            if self._model is None:
                raise Exception("not part of a model")
            if not hasattr(self._model, "solution") or \
                    self._model.solution is None or \
                    self._model.solution.status == "NA":
                raise Exception("model has not been solved")
            if self._model.solution.status != "optimal":
                raise Exception("model solution was not optimal")
            raise e  # Not sure what the exact problem was

    def shadow_price(self):
        """The shadow price for the metabolite in the most recent solution

        Shadow prices are computed from the dual values of the bounds in
        the solution.

        """
        try:
            return self._model.solution.shadow_prices[self.id]
        except Exception as e:
            if self._model is None:
                raise Exception("not part of a model")
            if not hasattr(self._model, "solution") or \
                    self._model.solution is None or \
                    self._model.solution.status == "NA":
                raise Exception("model has not been solved")
            if self._model.solution.status != "optimal":
                raise Exception("model solution was not optimal")
            raise e

    def remove_from_model(self, method='subtractive', **kwargs):
        """Removes the association from self.model

        Parameters
        ----------
        method : 'subtractive' or 'destructive'.
            If 'subtractive' then the metabolite is removed from all
            associated reactions.  If 'destructive' then all associated
            reactions are removed from the Model.
        """
        # why is model being taken in as a parameter? This plays
        # back to the question of allowing a Metabolite to be associated
        # with multiple Models
        if "model" in kwargs:
            warn("model argument deprecated")
        model = self._model
        self._model.metabolites.remove(self)
        self._model = None
        if method.lower() == 'subtractive':
            for the_reaction in list(self._reaction):
                the_coefficient = the_reaction._metabolites[self]
                the_reaction.subtract_metabolites({self: the_coefficient})
        elif method.lower() == 'destructive':
            for x in self._reaction:
                x.remove_from_model()
        else:
            raise Exception(method + " is not 'subtractive' or 'destructive'")
        model.solver.remove(model.solver.constraints[self.id])

    def summary(self, threshold=0.01, fva=False, floatfmt='.3g', **kwargs):
        """Print a summary of the reactions which produce and consume this
        metabolite. This method requires the model for which this metabolite is
        a part to be solved.

        Parameters
        ----------
        threshold : float
            a value below which to ignore reaction fluxes

        fva : float (0->1), or None
            Whether or not to include flux variability analysis in the output.
            If given, fva should be a float between 0 and 1, representing the
            fraction of the optimum objective to be searched.

        floatfmt : string
            format method for floats, passed to tabulate. Default is '.3g'.
        """
        try:
            from cobra.flux_analysis.summary import metabolite_summary
            return metabolite_summary(self, threshold=threshold, fva=fva,
                                      floatfmt=floatfmt, **kwargs)
        except ImportError:
            warn('Summary methods require pandas/tabulate')
