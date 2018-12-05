# -*- coding: utf-8 -*-

from __future__ import absolute_import, print_function

import hashlib
import re
from collections import defaultdict
from copy import copy, deepcopy
from functools import partial
from math import isinf
from operator import attrgetter
from warnings import warn

from future.utils import raise_from, raise_with_traceback
from six import iteritems, iterkeys, string_types

from cobra.core.configuration import Configuration
from cobra.core.gene import Gene, ast2str, eval_gpr, parse_gpr
from cobra.core.metabolite import Metabolite
from cobra.core.object import Object
from cobra.exceptions import OptimizationError
from cobra.util.context import get_context, resettable
from cobra.util.solver import (
    check_solver_status, linear_reaction_coefficients, set_objective)
from cobra.util.util import format_long_string


CONFIGURATION = Configuration()

# precompiled regular expressions
# Matches and/or in a gene reaction rule
and_or_search = re.compile(r'\(| and| or|\+|\)', re.IGNORECASE)
uppercase_AND = re.compile(r'\bAND\b')
uppercase_OR = re.compile(r'\bOR\b')
gpr_clean = re.compile(' {2,}')
# This regular expression finds any single letter compartment enclosed in
# square brackets at the beginning of the string. For example [c] : foo --> bar
compartment_finder = re.compile("^\s*(\[[A-Za-z]\])\s*:*")
# Regular expressions to match the arrows
_reversible_arrow_finder = re.compile("<(-+|=+)>")
_forward_arrow_finder = re.compile("(-+|=+)>")
_reverse_arrow_finder = re.compile("<(-+|=+)")


class Reaction(Object):
    """Reaction is a class for holding information regarding
    a biochemical reaction in a cobra.Model object.

    Parameters
    ----------
    id : string
        The identifier to associate with this reaction
    name : string
        A human readable name for the reaction
    subsystem : string
        Subsystem where the reaction is meant to occur
    lower_bound : float
        The lower flux bound
    upper_bound : float
        The upper flux bound
    """

    def __init__(self, id=None, name='', subsystem='', lower_bound=0.0,
                 upper_bound=None, objective_coefficient=0.0):
        Object.__init__(self, id, name)
        self._gene_reaction_rule = ''
        self.subsystem = subsystem

        # The cobra.Genes that are used to catalyze the reaction
        self._genes = set()

        # A dictionary of metabolites and their stoichiometric coefficients in
        # this reaction.
        self._metabolites = {}

        # The set of compartments that partaking metabolites are in.
        self._compartments = None

        # self.model is None or refers to the cobra.Model that
        # contains self
        self._model = None

        if objective_coefficient != 0:
            raise NotImplementedError('setting objective coefficient when '
                                      'creating reaction is no longer '
                                      'supported. Use the model.objective '
                                      'setter')

        # from cameo ...
        self._lower_bound = lower_bound if lower_bound is not None else \
            CONFIGURATION.lower_bound
        self._upper_bound = upper_bound if upper_bound is not None else \
            CONFIGURATION.upper_bound

    def _set_id_with_model(self, value):
        if value in self.model.reactions:
            raise ValueError("The model already contains a reaction with"
                             " the id:", value)
        forward_variable = self.forward_variable
        reverse_variable = self.reverse_variable
        self._id = value
        self.model.reactions._generate_index()
        forward_variable.name = self.id
        reverse_variable.name = self.reverse_id

    @property
    def reverse_id(self):
        """Generate the id of reverse_variable from the reaction's id."""
        return '_'.join((self.id, 'reverse',
                         hashlib.md5(
                             self.id.encode('utf-8')).hexdigest()[0:5]))

    @property
    def flux_expression(self):
        """Forward flux expression

        Returns
        -------
        sympy expression
            The expression representing the the forward flux (if associated
            with model), otherwise None. Representing the net flux if
            model.reversible_encoding == 'unsplit' or None if reaction is
            not associated with a model """
        if self.model is not None:
            return 1. * self.forward_variable - 1. * self.reverse_variable
        else:
            return None

    @property
    def forward_variable(self):
        """An optlang variable representing the forward flux

        Returns
        -------
        optlang.interface.Variable
            An optlang variable for the forward flux or None if reaction is
            not associated with a model.
        """
        if self.model is not None:
            return self.model.variables[self.id]
        else:
            return None

    @property
    def reverse_variable(self):
        """An optlang variable representing the reverse flux

        Returns
        -------
        optlang.interface.Variable
            An optlang variable for the reverse flux or None if reaction is
            not associated with a model.
        """

        if self.model is not None:
            return self.model.variables[self.reverse_id]
        else:
            return None

    @property
    def objective_coefficient(self):
        """ Get the coefficient for this reaction in a linear
        objective (float)

        Assuming that the objective of the associated model is summation of
        fluxes from a set of reactions, the coefficient for each reaction
        can be obtained individually using this property. A more general way
        is to use the `model.objective` property directly.
        """
        return linear_reaction_coefficients(self.model, [self]).get(self, 0)

    @objective_coefficient.setter
    def objective_coefficient(self, value):
        if self.model is None:
            raise AttributeError('cannot assign objective to a missing model')
        if self.flux_expression is not None:
            set_objective(self.model, {self: value}, additive=True)

    def __copy__(self):
        cop = copy(super(Reaction, self))
        return cop

    def __deepcopy__(self, memo):
        cop = deepcopy(super(Reaction, self), memo)
        return cop

    @staticmethod
    def _check_bounds(lb, ub):
        if lb > ub:
            raise ValueError(
                "The lower bound must be less than or equal to the upper "
                "bound ({} <= {}).".format(lb, ub))

    def update_variable_bounds(self):
        if self.model is None:
            return
        # We know that `lb <= ub`.
        if self._lower_bound > 0:
            self.forward_variable.set_bounds(
                lb=None if isinf(self._lower_bound) else self._lower_bound,
                ub=None if isinf(self._upper_bound) else self._upper_bound
            )
            self.reverse_variable.set_bounds(lb=0, ub=0)
        elif self._upper_bound < 0:
            self.forward_variable.set_bounds(lb=0, ub=0)
            self.reverse_variable.set_bounds(
                lb=None if isinf(self._upper_bound) else -self._upper_bound,
                ub=None if isinf(self._lower_bound) else -self._lower_bound
            )
        else:
            self.forward_variable.set_bounds(
                lb=0,
                ub=None if isinf(self._upper_bound) else self._upper_bound
            )
            self.reverse_variable.set_bounds(
                lb=0,
                ub=None if isinf(self._lower_bound) else -self._lower_bound
            )

    @property
    def lower_bound(self):
        """Get or set the lower bound

        Setting the lower bound (float) will also adjust the associated optlang
        variables associated with the reaction. Infeasible combinations,
        such as a lower bound higher than the current upper bound will
        update the other bound.

        When using a `HistoryManager` context, this attribute can be set
        temporarily, reversed when the exiting the context.
        """
        return self._lower_bound

    @lower_bound.setter
    @resettable
    def lower_bound(self, value):
        if self._upper_bound < value:
            warn("You are constraining the reaction '{}' to a fixed flux "
                 "value of {}. Did you intend to do this? We are planning to "
                 "remove this behavior in a future release. Please let us "
                 "know your opinion at "
                 "https://github.com/opencobra/cobrapy/issues/793."
                 "".format(self.id, value), DeprecationWarning)
            self.upper_bound = value
        # Validate bounds before setting them.
        self._check_bounds(value, self._upper_bound)
        self._lower_bound = value
        self.update_variable_bounds()

    @property
    def upper_bound(self):
        """Get or set the upper bound

        Setting the upper bound (float) will also adjust the associated optlang
        variables associated with the reaction. Infeasible combinations,
        such as a upper bound lower than the current lower bound will
        update the other bound.

        When using a `HistoryManager` context, this attribute can be set
        temporarily, reversed when the exiting the context.
        """
        return self._upper_bound

    @upper_bound.setter
    @resettable
    def upper_bound(self, value):
        if self._lower_bound > value:
            warn("You are constraining the reaction '{}' to a fixed flux "
                 "value of {}. Did you intend to do this? We are planning to "
                 "remove this behavior in a future release. Please let us "
                 "know your opinion at "
                 "https://github.com/opencobra/cobrapy/issues/793."
                 "".format(self.id, value), DeprecationWarning)
            self.lower_bound = value
        # Validate bounds before setting them.
        self._check_bounds(self._lower_bound, value)
        self._upper_bound = value
        self.update_variable_bounds()

    @property
    def bounds(self):
        """ Get or set the bounds directly from a tuple

        Convenience method for setting upper and lower bounds in one line
        using a tuple of lower and upper bound. Invalid bounds will raise an
        AssertionError.

        When using a `HistoryManager` context, this attribute can be set
        temporarily, reversed when the exiting the context.
        """
        return self.lower_bound, self.upper_bound

    @bounds.setter
    @resettable
    def bounds(self, value):
        lower, upper = value
        # Validate bounds before setting them.
        self._check_bounds(lower, upper)
        self._lower_bound = lower
        self._upper_bound = upper
        self.update_variable_bounds()

    @property
    def flux(self):
        """
        The flux value in the most recent solution.

        Flux is the primal value of the corresponding variable in the model.

        Warnings
        --------
        * Accessing reaction fluxes through a `Solution` object is the safer,
          preferred, and only guaranteed to be correct way. You can see how to
          do so easily in the examples.
        * Reaction flux is retrieved from the currently defined
          `self._model.solver`. The solver status is checked but there are no
          guarantees that the current solver state is the one you are looking
          for.
        * If you modify the underlying model after an optimization, you will
          retrieve the old optimization values.

        Raises
        ------
        RuntimeError
            If the underlying model was never optimized beforehand or the
            reaction is not part of a model.
        OptimizationError
            If the solver status is anything other than 'optimal'.
        AssertionError
            If the flux value is not within the bounds.

        Examples
        --------
        >>> import cobra.test
        >>> model = cobra.test.create_test_model("textbook")
        >>> solution = model.optimize()
        >>> model.reactions.PFK.flux
        7.477381962160283
        >>> solution.fluxes.PFK
        7.4773819621602833
        """
        try:
            check_solver_status(self._model.solver.status)
            return self.forward_variable.primal - self.reverse_variable.primal
        except AttributeError:
            raise RuntimeError(
                "reaction '{}' is not part of a model".format(self.id))
        # Due to below all-catch, which sucks, need to reraise these.
        except (RuntimeError, OptimizationError) as err:
            raise_with_traceback(err)
        # Would love to catch CplexSolverError and GurobiError here.
        except Exception as err:
            raise_from(OptimizationError(
                "Likely no solution exists. Original solver message: {}."
                "".format(str(err))), err)

    @property
    def reduced_cost(self):
        """
        The reduced cost in the most recent solution.

        Reduced cost is the dual value of the corresponding variable in the
        model.

        Warnings
        --------
        * Accessing reduced costs through a `Solution` object is the safer,
          preferred, and only guaranteed to be correct way. You can see how to
          do so easily in the examples.
        * Reduced cost is retrieved from the currently defined
          `self._model.solver`. The solver status is checked but there are no
          guarantees that the current solver state is the one you are looking
          for.
        * If you modify the underlying model after an optimization, you will
          retrieve the old optimization values.

        Raises
        ------
        RuntimeError
            If the underlying model was never optimized beforehand or the
            reaction is not part of a model.
        OptimizationError
            If the solver status is anything other than 'optimal'.

        Examples
        --------
        >>> import cobra.test
        >>> model = cobra.test.create_test_model("textbook")
        >>> solution = model.optimize()
        >>> model.reactions.PFK.reduced_cost
        -8.673617379884035e-18
        >>> solution.reduced_costs.PFK
        -8.6736173798840355e-18
        """
        try:
            check_solver_status(self._model.solver.status)
            return self.forward_variable.dual - self.reverse_variable.dual
        except AttributeError:
            raise RuntimeError(
                "reaction '{}' is not part of a model".format(self.id))
        # Due to below all-catch, which sucks, need to reraise these.
        except (RuntimeError, OptimizationError) as err:
            raise_with_traceback(err)
        # Would love to catch CplexSolverError and GurobiError here.
        except Exception as err:
            raise_from(OptimizationError(
                "Likely no solution exists. Original solver message: {}."
                "".format(str(err))), err)

    # read-only
    @property
    def metabolites(self):
        return self._metabolites.copy()

    @property
    def genes(self):
        return frozenset(self._genes)

    @property
    def gene_reaction_rule(self):
        return self._gene_reaction_rule

    @gene_reaction_rule.setter
    def gene_reaction_rule(self, new_rule):

        # TODO: Do this :)
        if get_context(self):
            warn("Context management not implemented for "
                 "gene reaction rules")

        self._gene_reaction_rule = new_rule.strip()
        try:
            _, gene_names = parse_gpr(self._gene_reaction_rule)
        except (SyntaxError, TypeError) as e:
            if "AND" in new_rule or "OR" in new_rule:
                warn("uppercase AND/OR found in rule '%s' for '%s'" %
                     (new_rule, repr(self)))
                new_rule = uppercase_AND.sub("and", new_rule)
                new_rule = uppercase_OR.sub("or", new_rule)
                self.gene_reaction_rule = new_rule
                return
            warn("malformed gene_reaction_rule '%s' for %s" %
                 (new_rule, repr(self)))
            tmp_str = and_or_search.sub('', self._gene_reaction_rule)
            gene_names = set((gpr_clean.sub(' ', tmp_str).split(' ')))
        if '' in gene_names:
            gene_names.remove('')
        old_genes = self._genes
        if self._model is None:
            self._genes = {Gene(i) for i in gene_names}
        else:
            model_genes = self._model.genes
            self._genes = set()
            for id in gene_names:
                if model_genes.has_id(id):
                    self._genes.add(model_genes.get_by_id(id))
                else:
                    new_gene = Gene(id)
                    new_gene._model = self._model
                    self._genes.add(new_gene)
                    model_genes.append(new_gene)

        # Make the genes aware that it is involved in this reaction
        for g in self._genes:
            g._reaction.add(self)

        # make the old genes aware they are no longer involved in this reaction
        for g in old_genes:
            if g not in self._genes:  # if an old gene is not a new gene
                try:
                    g._reaction.remove(self)
                except KeyError:
                    warn("could not remove old gene %s from reaction %s" %
                         (g.id, self.id))

    @property
    def gene_name_reaction_rule(self):
        """Display gene_reaction_rule with names intead.

        Do NOT use this string for computation. It is intended to give a
        representation of the rule using more familiar gene names instead of
        the often cryptic ids.

        """
        names = {i.id: i.name for i in self._genes}
        ast = parse_gpr(self._gene_reaction_rule)[0]
        return ast2str(ast, names=names)

    @property
    def functional(self):
        """All required enzymes for reaction are functional.

        Returns
        -------
        bool
            True if the gene-protein-reaction (GPR) rule is fulfilled for
            this reaction, or if reaction is not associated to a model,
            otherwise False.
        """
        if self._model:
            tree, _ = parse_gpr(self.gene_reaction_rule)
            return eval_gpr(tree, {gene.id for gene in self.genes if
                                   not gene.functional})
        return True

    @property
    def x(self):
        """The flux through the reaction in the most recent solution.

        Flux values are computed from the primal values of the variables in
        the solution.
        """
        warn("Please use reaction.flux instead.", DeprecationWarning)
        return self.flux

    @property
    def y(self):
        """The reduced cost of the reaction in the most recent solution.

        Reduced costs are computed from the dual values of the variables in
        the solution.
        """
        warn("Please use reaction.reduced_cost instead.", DeprecationWarning)
        return self.reduced_cost

    @property
    def reversibility(self):
        """Whether the reaction can proceed in both directions (reversible)

        This is computed from the current upper and lower bounds.

        """
        return self._lower_bound < 0 < self._upper_bound

    @reversibility.setter
    def reversibility(self, value):
        warn("Setting reaction reversibility is ignored")

    @property
    def boundary(self):
        """Whether or not this reaction is an exchange reaction.

        Returns `True` if the reaction has either no products or reactants.
        """
        return (len(self.metabolites) == 1 and
                not (self.reactants and self.products))

    @property
    def model(self):
        """returns the model the reaction is a part of"""
        return self._model

    def _update_awareness(self):
        """Make sure all metabolites and genes that are associated with
        this reaction are aware of it.

        """
        for x in self._metabolites:
            x._reaction.add(self)
        for x in self._genes:
            x._reaction.add(self)

    def remove_from_model(self, remove_orphans=False):
        """Removes the reaction from a model.

        This removes all associations between a reaction the associated
        model, metabolites and genes.

        The change is reverted upon exit when using the model as a context.

        Parameters
        ----------
        remove_orphans : bool
            Remove orphaned genes and metabolites from the model as well

        """
        self._model.remove_reactions([self], remove_orphans=remove_orphans)

    def delete(self, remove_orphans=False):
        """Removes the reaction from a model.

        This removes all associations between a reaction the associated
        model, metabolites and genes.

        The change is reverted upon exit when using the model as a context.

        Deprecated, use `reaction.remove_from_model` instead.

        Parameters
        ----------
        remove_orphans : bool
            Remove orphaned genes and metabolites from the model as well

        """
        warn("delete is deprecated. Use reaction.remove_from_model instead",
             DeprecationWarning)
        self.remove_from_model(remove_orphans=remove_orphans)

    def __setstate__(self, state):
        """Probably not necessary to set _model as the cobra.Model that
        contains self sets the _model attribute for all metabolites and genes
        in the reaction.

        However, to increase performance speed we do want to let the metabolite
        and gene know that they are employed in this reaction

        """
        # These are necessary for old pickles which store attributes
        # which have since been superceded by properties.
        if "reaction" in state:
            state.pop("reaction")
        if "gene_reaction_rule" in state:
            state["_gene_reaction_rule"] = state.pop("gene_reaction_rule")
        if "lower_bound" in state:
            state['_lower_bound'] = state.pop('lower_bound')
        if "upper_bound" in state:
            state['_upper_bound'] = state.pop('upper_bound')

        self.__dict__.update(state)
        for x in state['_metabolites']:
            setattr(x, '_model', self._model)
            x._reaction.add(self)
        for x in state['_genes']:
            setattr(x, '_model', self._model)
            x._reaction.add(self)

    def copy(self):
        """Copy a reaction

        The referenced metabolites and genes are also copied.

        """
        # no references to model when copying
        model = self._model
        self._model = None
        for i in self._metabolites:
            i._model = None
        for i in self._genes:
            i._model = None
        # now we can copy
        new_reaction = deepcopy(self)
        # restore the references
        self._model = model
        for i in self._metabolites:
            i._model = model
        for i in self._genes:
            i._model = model
        return new_reaction

    def __add__(self, other):
        """Add two reactions

        The stoichiometry will be the combined stoichiometry of the two
        reactions, and the gene reaction rule will be both rules combined by an
        and. All other attributes (i.e. reaction bounds) will match those of
        the first reaction

        """
        new_reaction = self.copy()
        if other == 0:
            return new_reaction
        else:
            new_reaction += other

        return new_reaction

    __radd__ = __add__

    def __iadd__(self, other):

        self.add_metabolites(other._metabolites, combine=True)
        gpr1 = self.gene_reaction_rule.strip()
        gpr2 = other.gene_reaction_rule.strip()
        if gpr1 != '' and gpr2 != '':
            self.gene_reaction_rule = "(%s) and (%s)" % \
                                      (self.gene_reaction_rule,
                                       other.gene_reaction_rule)
        elif gpr1 != '' and gpr2 == '':
            self.gene_reaction_rule = gpr1
        elif gpr1 == '' and gpr2 != '':
            self.gene_reaction_rule = gpr2
        return self

    def __sub__(self, other):
        new = self.copy()
        new -= other
        return new

    def __isub__(self, other):
        self.subtract_metabolites(other._metabolites, combine=True)
        return self

    def __imul__(self, coefficient):
        """Scale coefficients in a reaction by a given value

        E.g. A -> B becomes 2A -> 2B.

        If coefficient is less than zero, the reaction is reversed and the
        bounds are swapped.
        """
        self._metabolites = {
            met: value * coefficient
            for met, value in iteritems(self._metabolites)}

        if coefficient < 0:
            self.bounds = (-self.upper_bound, -self.lower_bound)

        if self._model:
            self._model._populate_solver([self])

        context = get_context(self)
        if context:
            context(partial(self._model._populate_solver, [self]))
            context(partial(self.__imul__, 1./coefficient))

        return self

    def __mul__(self, coefficient):
        new = self.copy()
        new *= coefficient
        return new

    @property
    def reactants(self):
        """Return a list of reactants for the reaction."""
        return [k for k, v in iteritems(self._metabolites) if v < 0]

    @property
    def products(self):
        """Return a list of products for the reaction"""
        return [k for k, v in iteritems(self._metabolites) if v >= 0]

    def get_coefficient(self, metabolite_id):
        """
        Return the stoichiometric coefficient of a metabolite.

        Parameters
        ----------
        metabolite_id : str or cobra.Metabolite

        """
        if isinstance(metabolite_id, Metabolite):
            return self._metabolites[metabolite_id]

        _id_to_metabolites = {m.id: m for m in self._metabolites}
        return self._metabolites[_id_to_metabolites[metabolite_id]]

    def get_coefficients(self, metabolite_ids):
        """
        Return the stoichiometric coefficients for a list of metabolites.

        Parameters
        ----------
        metabolite_ids : iterable
            Containing ``str`` or ``cobra.Metabolite``s.

        """
        return map(self.get_coefficient, metabolite_ids)

    def add_metabolites(self, metabolites_to_add, combine=True,
                        reversibly=True):
        """Add metabolites and stoichiometric coefficients to the reaction.
        If the final coefficient for a metabolite is 0 then it is removed
        from the reaction.

        The change is reverted upon exit when using the model as a context.

        Parameters
        ----------
        metabolites_to_add : dict
            Dictionary with metabolite objects or metabolite identifiers as
            keys and coefficients as values. If keys are strings (name of a
            metabolite) the reaction must already be part of a model and a
            metabolite with the given name must exist in the model.

        combine : bool
            Describes behavior a metabolite already exists in the reaction.
            True causes the coefficients to be added.
            False causes the coefficient to be replaced.

        reversibly : bool
            Whether to add the change to the context to make the change
            reversibly or not (primarily intended for internal use).

        """
        old_coefficients = self.metabolites
        new_metabolites = []
        _id_to_metabolites = dict([(x.id, x) for x in self._metabolites])

        for metabolite, coefficient in iteritems(metabolites_to_add):

            # Make sure metabolites being added belong to the same model, or
            # else copy them.
            if isinstance(metabolite, Metabolite):
                if ((metabolite.model is not None) and
                        (metabolite.model is not self._model)):
                    metabolite = metabolite.copy()

            met_id = str(metabolite)
            # If a metabolite already exists in the reaction then
            # just add them.
            if met_id in _id_to_metabolites:
                reaction_metabolite = _id_to_metabolites[met_id]
                if combine:
                    self._metabolites[reaction_metabolite] += coefficient
                else:
                    self._metabolites[reaction_metabolite] = coefficient
            else:
                # If the reaction is in a model, ensure we aren't using
                # a duplicate metabolite.
                if self._model:
                    try:
                        metabolite = \
                            self._model.metabolites.get_by_id(met_id)
                    except KeyError as e:
                        if isinstance(metabolite, Metabolite):
                            new_metabolites.append(metabolite)
                        else:
                            # do we want to handle creation here?
                            raise e
                elif isinstance(metabolite, string_types):
                    # if we want to handle creation, this should be changed
                    raise ValueError("Reaction '%s' does not belong to a "
                                     "model. Either add the reaction to a "
                                     "model or use Metabolite objects instead "
                                     "of strings as keys."
                                     % self.id)
                self._metabolites[metabolite] = coefficient
                # make the metabolite aware that it is involved in this
                # reaction
                metabolite._reaction.add(self)

        # from cameo ...
        model = self.model
        if model is not None:
            model.add_metabolites(new_metabolites)

            for metabolite, coefficient in self._metabolites.items():
                model.constraints[
                    metabolite.id].set_linear_coefficients(
                    {self.forward_variable: coefficient,
                     self.reverse_variable: -coefficient
                     })

        for metabolite, the_coefficient in list(self._metabolites.items()):
            if the_coefficient == 0:
                # make the metabolite aware that it no longer participates
                # in this reaction
                metabolite._reaction.remove(self)
                self._metabolites.pop(metabolite)

        context = get_context(self)
        if context and reversibly:
            if combine:
                # Just subtract the metabolites that were added
                context(partial(
                    self.subtract_metabolites, metabolites_to_add,
                    combine=True, reversibly=False))
            else:
                # Reset them with add_metabolites
                mets_to_reset = {
                    key: old_coefficients[model.metabolites.get_by_any(key)[0]]
                    for key in iterkeys(metabolites_to_add)}

                context(partial(
                    self.add_metabolites, mets_to_reset,
                    combine=False, reversibly=False))

    def subtract_metabolites(self, metabolites, combine=True, reversibly=True):
        """Subtract metabolites from a reaction.

        That means add the metabolites with -1*coefficient. If the final
        coefficient for a metabolite is 0 then the metabolite is removed from
        the reaction.

        Notes
        -----
        * A final coefficient < 0 implies a reactant.
        * The change is reverted upon exit when using the model as a context.

        Parameters
        ----------
        metabolites : dict
            Dictionary where the keys are of class Metabolite and the values
            are the coefficients. These metabolites will be added to the
            reaction.

        combine : bool
            Describes behavior a metabolite already exists in the reaction.
            True causes the coefficients to be added.
            False causes the coefficient to be replaced.

        reversibly : bool
            Whether to add the change to the context to make the change
            reversibly or not (primarily intended for internal use).

        """
        self.add_metabolites({
            k: -v for k, v in iteritems(metabolites)},
            combine=combine, reversibly=reversibly)

    @property
    def reaction(self):
        """Human readable reaction string"""
        return self.build_reaction_string()

    @reaction.setter
    def reaction(self, value):
        return self.build_reaction_from_string(value)

    def build_reaction_string(self, use_metabolite_names=False):
        """Generate a human readable reaction string"""

        def format(number):
            return "" if number == 1 else str(number).rstrip(".") + " "

        id_type = 'id'
        if use_metabolite_names:
            id_type = 'name'
        reactant_bits = []
        product_bits = []
        for met in sorted(self._metabolites, key=attrgetter("id")):
            coefficient = self._metabolites[met]
            name = str(getattr(met, id_type))
            if coefficient >= 0:
                product_bits.append(format(coefficient) + name)
            else:
                reactant_bits.append(format(abs(coefficient)) + name)

        reaction_string = ' + '.join(reactant_bits)
        if not self.reversibility:
            if self.lower_bound < 0 and self.upper_bound <= 0:
                reaction_string += ' <-- '
            else:
                reaction_string += ' --> '
        else:
            reaction_string += ' <=> '
        reaction_string += ' + '.join(product_bits)
        return reaction_string

    def check_mass_balance(self):
        """Compute mass and charge balance for the reaction

        returns a dict of {element: amount} for unbalanced elements.
        "charge" is treated as an element in this dict
        This should be empty for balanced reactions.
        """
        reaction_element_dict = defaultdict(int)
        for metabolite, coefficient in iteritems(self._metabolites):
            if metabolite.charge is not None:
                reaction_element_dict["charge"] += \
                    coefficient * metabolite.charge
            if metabolite.elements is None:
                raise ValueError("No elements found in metabolite %s"
                                 % metabolite.id)
            for element, amount in iteritems(metabolite.elements):
                reaction_element_dict[element] += coefficient * amount
        # filter out 0 values
        return {k: v for k, v in iteritems(reaction_element_dict) if v != 0}

    @property
    def compartments(self):
        """lists compartments the metabolites are in"""
        if self._compartments is None:
            self._compartments = {met.compartment for met in self._metabolites
                                  if met.compartment is not None}
        return self._compartments

    def get_compartments(self):
        """lists compartments the metabolites are in"""
        warn('use Reaction.compartments instead', DeprecationWarning)
        return list(self.compartments)

    def _associate_gene(self, cobra_gene):
        """Associates a cobra.Gene object with a cobra.Reaction.

        Parameters
        ----------
        cobra_gene : cobra.core.Gene.Gene

        """
        self._genes.add(cobra_gene)
        cobra_gene._reaction.add(self)
        cobra_gene._model = self._model

    def _dissociate_gene(self, cobra_gene):
        """Dissociates a cobra.Gene object with a cobra.Reaction.

        Parameters
        ----------
        cobra_gene : cobra.core.Gene.Gene

        """
        self._genes.discard(cobra_gene)
        cobra_gene._reaction.discard(self)

    def knock_out(self):
        """Knockout reaction by setting its bounds to zero."""
        self.bounds = (0, 0)

    def build_reaction_from_string(self, reaction_str, verbose=True,
                                   fwd_arrow=None, rev_arrow=None,
                                   reversible_arrow=None, term_split="+"):
        """Builds reaction from reaction equation reaction_str using parser

        Takes a string and using the specifications supplied in the optional
        arguments infers a set of metabolites, metabolite compartments and
        stoichiometries for the reaction.  It also infers the reversibility
        of the reaction from the reaction arrow.

        Changes to the associated model are reverted upon exit when using
        the model as a context.

        Parameters
        ----------
        reaction_str : string
            a string containing a reaction formula (equation)
        verbose: bool
            setting verbosity of function
        fwd_arrow : re.compile
            for forward irreversible reaction arrows
        rev_arrow : re.compile
            for backward irreversible reaction arrows
        reversible_arrow : re.compile
            for reversible reaction arrows
        term_split : string
            dividing individual metabolite entries

        """
        # set the arrows
        forward_arrow_finder = _forward_arrow_finder if fwd_arrow is None \
            else re.compile(re.escape(fwd_arrow))
        reverse_arrow_finder = _reverse_arrow_finder if rev_arrow is None \
            else re.compile(re.escape(rev_arrow))
        reversible_arrow_finder = _reversible_arrow_finder \
            if reversible_arrow is None \
            else re.compile(re.escape(reversible_arrow))
        if self._model is None:
            warn("no model found")
            model = None
        else:
            model = self._model
        found_compartments = compartment_finder.findall(reaction_str)
        if len(found_compartments) == 1:
            compartment = found_compartments[0]
            reaction_str = compartment_finder.sub("", reaction_str)
        else:
            compartment = ""

        # reversible case
        arrow_match = reversible_arrow_finder.search(reaction_str)
        if arrow_match is not None:
            self.lower_bound = -1000
            self.upper_bound = 1000
        else:  # irreversible
            # try forward
            arrow_match = forward_arrow_finder.search(reaction_str)
            if arrow_match is not None:
                self.upper_bound = 1000
                self.lower_bound = 0
            else:
                # must be reverse
                arrow_match = reverse_arrow_finder.search(reaction_str)
                if arrow_match is None:
                    raise ValueError("no suitable arrow found in '%s'" %
                                     reaction_str)
                else:
                    self.upper_bound = 0
                    self.lower_bound = -1000
        reactant_str = reaction_str[:arrow_match.start()].strip()
        product_str = reaction_str[arrow_match.end():].strip()

        self.subtract_metabolites(self.metabolites, combine=True)

        for substr, factor in ((reactant_str, -1), (product_str, 1)):
            if len(substr) == 0:
                continue
            for term in substr.split(term_split):
                term = term.strip()
                if term.lower() == "nothing":
                    continue
                if " " in term:
                    num_str, met_id = term.split()
                    num = float(num_str.lstrip("(").rstrip(")")) * factor
                else:
                    met_id = term
                    num = factor
                met_id += compartment
                try:
                    met = model.metabolites.get_by_id(met_id)
                except KeyError:
                    if verbose:
                        print("unknown metabolite '%s' created" % met_id)
                    met = Metabolite(met_id)
                self.add_metabolites({met: num})

    def __str__(self):
        return "{id}: {stoichiometry}".format(
            id=self.id, stoichiometry=self.build_reaction_string())

    def _repr_html_(self):
        return """
        <table>
            <tr>
                <td><strong>Reaction identifier</strong></td><td>{id}</td>
            </tr><tr>
                <td><strong>Name</strong></td><td>{name}</td>
            </tr><tr>
                <td><strong>Memory address</strong></td>
                <td>{address}</td>
            </tr><tr>
                <td><strong>Stoichiometry</strong></td>
                <td>
                    <p style='text-align:right'>{stoich_id}</p>
                    <p style='text-align:right'>{stoich_name}</p>
                </td>
            </tr><tr>
                <td><strong>GPR</strong></td><td>{gpr}</td>
            </tr><tr>
                <td><strong>Lower bound</strong></td><td>{lb}</td>
            </tr><tr>
                <td><strong>Upper bound</strong></td><td>{ub}</td>
            </tr>
        </table>
        """.format(id=format_long_string(self.id, 100),
                   name=format_long_string(self.name, 100),
                   address='0x0%x' % id(self),
                   stoich_id=format_long_string(
                       self.build_reaction_string(), 200),
                   stoich_name=format_long_string(
                       self.build_reaction_string(True), 200),
                   gpr=format_long_string(self.gene_reaction_rule, 100),
                   lb=self.lower_bound, ub=self.upper_bound)
