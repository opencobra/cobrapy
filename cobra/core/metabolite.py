# -*- coding: utf-8 -*-

from __future__ import absolute_import

import re
from warnings import warn

from future.utils import raise_from, raise_with_traceback
from six import iteritems

from cobra.core.formula import elements_and_molecular_weights
from cobra.core.species import Species
from cobra.exceptions import OptimizationError
from cobra.util.solver import check_solver_status
from cobra.util.util import format_long_string


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

        self._bound = 0.

    def _set_id_with_model(self, value):
        if value in self.model.metabolites:
            raise ValueError("The model already contains a metabolite with "
                             "the id:", value)
        self.model.constraints[self.id].name = value
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
            return self.model.constraints[self.id]

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
        # commonly occurring characters in incorrectly constructed formulas
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
            warn("The element %s does not appear in the periodic table" % e)

    @property
    def y(self):
        """The shadow price for the metabolite in the most recent solution

        Shadow prices are computed from the dual values of the bounds in
        the solution.

        """
        warn("Please use metabolite.shadow_price instead.", DeprecationWarning)
        return self.shadow_price

    @property
    def shadow_price(self):
        """
        The shadow price in the most recent solution.

        Shadow price is the dual value of the corresponding constraint in the
        model.

        Warnings
        --------
        * Accessing shadow prices through a `Solution` object is the safer,
          preferred, and only guaranteed to be correct way. You can see how to
          do so easily in the examples.
        * Shadow price is retrieved from the currently defined
          `self._model.solver`. The solver status is checked but there are no
          guarantees that the current solver state is the one you are looking
          for.
        * If you modify the underlying model after an optimization, you will
          retrieve the old optimization values.

        Raises
        ------
        RuntimeError
            If the underlying model was never optimized beforehand or the
            metabolite is not part of a model.
        OptimizationError
            If the solver status is anything other than 'optimal'.

        Examples
        --------
        >>> import cobra
        >>> import cobra.test
        >>> model = cobra.test.create_test_model("textbook")
        >>> solution = model.optimize()
        >>> model.metabolites.glc__D_e.shadow_price
        -0.09166474637510488
        >>> solution.shadow_prices.glc__D_e
        -0.091664746375104883
        """
        try:
            check_solver_status(self._model.solver.status)
            return self._model.constraints[self.id].dual
        except AttributeError:
            raise RuntimeError(
                "metabolite '{}' is not part of a model".format(self.id))
        # Due to below all-catch, which sucks, need to reraise these.
        except (RuntimeError, OptimizationError) as err:
            raise_with_traceback(err)
        # Would love to catch CplexSolverError and GurobiError here.
        except Exception as err:
            raise_from(OptimizationError(
                "Likely no solution exists. Original solver message: {}."
                "".format(str(err))), err)

    def remove_from_model(self, destructive=False):
        """Removes the association from self.model

        The change is reverted upon exit when using the model as a context.

        Parameters
        ----------
        destructive : bool
            If False then the metabolite is removed from all
            associated reactions.  If True then all associated
            reactions are removed from the Model.
        """
        self._model.remove_metabolites(self, destructive)

    def summary(self, solution=None, threshold=0.01, fva=None, names=False,
                floatfmt='.3g'):
        """
        Print a summary of the production and consumption fluxes.

        This method requires the model for which this metabolite is a part
        to be solved.

        Parameters
        ----------
        solution : cobra.Solution, optional
            A previously solved model solution to use for generating the
            summary. If none provided (default), the summary method will
            resolve the model. Note that the solution object must match the
            model, i.e., changes to the model such as changed bounds,
            added or removed reactions are not taken into account by this
            method.
        threshold : float, optional
            Threshold below which fluxes are not reported.
        fva : pandas.DataFrame, float or None, optional
            Whether or not to include flux variability analysis in the output.
            If given, fva should either be a previous FVA solution matching
            the model or a float between 0 and 1 representing the
            fraction of the optimum objective to be searched.
        names : bool, optional
            Emit reaction and metabolite names rather than identifiers (default
            False).
        floatfmt : string, optional
            Format string for floats (default '.3g').

        """
        from cobra.flux_analysis.summary import metabolite_summary
        return metabolite_summary(self, solution=solution, threshold=threshold,
                                  fva=fva, names=names, floatfmt=floatfmt)

    def _repr_html_(self):
        return """
        <table>
            <tr>
                <td><strong>Metabolite identifier</strong></td><td>{id}</td>
            </tr><tr>
                <td><strong>Name</strong></td><td>{name}</td>
            </tr><tr>
                <td><strong>Memory address</strong></td>
                <td>{address}</td>
            </tr><tr>
                <td><strong>Formula</strong></td><td>{formula}</td>
            </tr><tr>
                <td><strong>Compartment</strong></td><td>{compartment}</td>
            </tr><tr>
                <td><strong>In {n_reactions} reaction(s)</strong></td><td>
                    {reactions}</td>
            </tr>
        </table>""".format(id=self.id, name=format_long_string(self.name),
                           formula=self.formula,
                           address='0x0%x' % id(self),
                           compartment=self.compartment,
                           n_reactions=len(self.reactions),
                           reactions=format_long_string(
                               ', '.join(r.id for r in self.reactions), 200))
