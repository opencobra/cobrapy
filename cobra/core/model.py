# -*- coding: utf-8 -*-

from __future__ import absolute_import

import types
from copy import copy, deepcopy
from functools import partial
from warnings import warn

import optlang
import six
import sympy
from six import iteritems, string_types
from sympy import S

from cobra.core.dictlist import DictList
from cobra.core.object import Object
from cobra.core.reaction import separate_forward_and_reverse_bounds
from cobra.core.solution import get_solution
from cobra.solvers import optimize
from cobra.util.context import HistoryManager, resettable, get_context
from cobra.util.solver import (
    SolverNotFound, get_solver_name, interface_to_str, set_objective, solvers,
    add_cons_vars_to_problem, remove_cons_vars_from_problem, choose_solver,
    check_solver_status)
from cobra.util.util import AutoVivification


class Model(Object):
    """Class representation for a cobra model

    Parameters
    ----------
    id_or_model : Model, string
        Either an existing Model object in which case a new model object is
        instantiated with the same properties as the original model,
        or a the identifier to associate with the model as a string.
    name : string
        Human readable name for the model

    Attributes
    ----------
    reactions : DictList
        A DictList where the key is the reaction identifier and the value a
        Reaction
    metabolites : DictList
        A DictList where the key is the metabolite identifier and the value a
        Metabolite
    genes : DictList
        A DictList where the key is the gene identifier and the value a
        Gene
    solution : Solution
        The last obtained solution from optimizing the model.
    """

    def __setstate__(self, state):
        """Make sure all cobra.Objects in the model point to the model.
        """
        self.__dict__.update(state)
        for y in ['reactions', 'genes', 'metabolites']:
            for x in getattr(self, y):
                x._model = self
                if y == 'reactions':
                    x._reset_var_cache()
        if not hasattr(self, "name"):
            self.name = None

    def __getstate__(self):
        """Get state for serialization.

        Ensures that the context stack is cleared prior to serialization,
        since partial functions cannot be pickled reliably.
        """
        odict = self.__dict__.copy()
        odict['_contexts'] = []
        return odict

    def __init__(self, id_or_model=None, name=None):
        if isinstance(id_or_model, Model):
            Object.__init__(self, name=name)
            self.__setstate__(id_or_model.__dict__)
            if not hasattr(self, "name"):
                self.name = None
            self._solver = id_or_model.solver
        else:
            Object.__init__(self, id_or_model, name=name)
            self._trimmed = False
            self._trimmed_genes = []
            self._trimmed_reactions = {}
            self.genes = DictList()
            self.reactions = DictList()  # A list of cobra.Reactions
            self.metabolites = DictList()  # A list of cobra.Metabolites
            # genes based on their ids {Gene.id: Gene}
            self.compartments = dict()
            self._contexts = []

            # from cameo ...

            # if not hasattr(self, '_solver'):  # backwards compatibility
            # with older cobrapy pickles?
            interface = solvers[get_solver_name()]
            self._solver = interface.Model()
            self._solver.objective = interface.Objective(S.Zero)
            self._populate_solver(self.reactions, self.metabolites)

    @property
    def solver(self):
        """Get or set the attached solver instance.

        The associated the solver object, which manages the interaction with
        the associated solver, e.g. glpk.

        This property is useful for accessing the optimization problem
        directly and to define additional non-metabolic constraints.

        Examples
        --------
        >>> import cobra.test
        >>> model = cobra.test.create_test_model("textbook")
        >>> new = model.problem.Constraint(model.objective.expression,
        >>> lb=0.99)
        >>> model.solver.add(new)
        """
        return self._solver

    @solver.setter
    @resettable
    def solver(self, value):
        not_valid_interface = SolverNotFound(
            '%s is not a valid solver interface. Pick from %s, or specify an '
            'optlang interface (e.g. optlang.glpk_interface).' % (
                value, list(solvers.keys())))
        if isinstance(value, six.string_types):
            try:
                interface = solvers[interface_to_str(value)]
            except KeyError:
                raise not_valid_interface
        elif isinstance(value, types.ModuleType) and hasattr(value, 'Model'):
            interface = value
        elif isinstance(value, optlang.interface.Model):
            interface = value.interface
        else:
            raise not_valid_interface

        # Do nothing if the solver did not change
        if self.problem == interface:
            return

        for reaction in self.reactions:
            reaction._reset_var_cache()
        self._solver = interface.Model.clone(self._solver)

    @property
    def description(self):
        warn("description deprecated", DeprecationWarning)
        return self.name if self.name is not None else ""

    @description.setter
    def description(self, value):
        self.name = value
        warn("description deprecated", DeprecationWarning)

    def get_metabolite_compartments(self):
        """Return all metabolites' compartments."""
        return {met.compartment for met in self.metabolites
                if met.compartment is not None}

    @property
    def medium(self):

        def is_active(reaction):
            """Determine if a boundary reaction permits flux towards creating
            metabolites
            """

            return ((bool(reaction.products) and (reaction.upper_bound > 0)) or
                    (bool(reaction.reactants) and (reaction.lower_bound < 0)))

        def get_active_bound(reaction):
            """For an active boundary reaction, return the relevant bound"""
            if reaction.reactants:
                return -reaction.lower_bound
            elif reaction.products:
                return reaction.upper_bound

        return {rxn.id: get_active_bound(rxn) for rxn in self.exchanges
                if is_active(rxn)}

    @medium.setter
    def medium(self, medium):
        """Get or set the constraints on the model exchanges.

        `model.medium` returns a dictionary of the bounds for each of the
        boundary reactions, in the form of `{rxn_id: bound}`, where `bound`
        specifies the absolute value of the bound in direction of metabolite
        creation (i.e., lower_bound for `met <--`, upper_bound for `met -->`)

        Parameters
        ----------
        medium: dictionary-like
            The medium to initialize. medium should be a dictionary defining
            `{rxn_id: bound}` pairs.

        """

        def set_active_bound(reaction, bound):
            if reaction.reactants:
                reaction.lower_bound = -bound
            elif reaction.products:
                reaction.upper_bound = bound

        # Set the given media bounds
        media_rxns = list()
        for rxn_id, bound in iteritems(medium):
            rxn = self.reactions.get_by_id(rxn_id)
            media_rxns.append(rxn)
            set_active_bound(rxn, bound)

        boundary_rxns = set(self.exchanges)
        media_rxns = set(media_rxns)

        # Turn off reactions not present in media
        for rxn in (boundary_rxns - media_rxns):
            set_active_bound(rxn, 0)

    def __add__(self, other_model):
        """Adds two models. +

        The issue of reactions being able to exists in multiple Models now
        arises, the same for metabolites and such.  This might be a little
        difficult as a reaction with the same name / id in two models might
        have different coefficients for their metabolites due to pH and whatnot
        making them different reactions.

        """
        new_model = self.copy()
        new_reactions = deepcopy(other_model.reactions)
        new_model.add_reactions(new_reactions)
        new_model.id = self.id + '_' + other_model.id
        return new_model

    def __iadd__(self, other_model):
        """Adds a Model to this model +=

        The issue of reactions being able to exists in multiple Models now
        arises, the same for metabolites and such.  This might be a little
        difficult as a reaction with the same name / id in two models might
        have different coefficients for their metabolites due to pH and whatnot
        making them different reactions.

        """
        new_reactions = deepcopy(other_model.reactions)
        self.add_reactions(new_reactions)
        self.id = self.id + '_' + other_model.id
        return self

    def copy(self):
        """Provides a partial 'deepcopy' of the Model.  All of the Metabolite,
        Gene, and Reaction objects are created anew but in a faster fashion
        than deepcopy
        """
        new = self.__class__()
        do_not_copy_by_ref = {"metabolites", "reactions", "genes", "notes",
                              "annotation"}
        for attr in self.__dict__:
            if attr not in do_not_copy_by_ref:
                new.__dict__[attr] = self.__dict__[attr]
        new.notes = deepcopy(self.notes)
        new.annotation = deepcopy(self.annotation)

        new.metabolites = DictList()
        do_not_copy_by_ref = {"_reaction", "_model"}
        for metabolite in self.metabolites:
            new_met = metabolite.__class__()
            for attr, value in iteritems(metabolite.__dict__):
                if attr not in do_not_copy_by_ref:
                    new_met.__dict__[attr] = copy(
                        value) if attr == "formula" else value
            new_met._model = new
            new.metabolites.append(new_met)

        new.genes = DictList()
        for gene in self.genes:
            new_gene = gene.__class__(None)
            for attr, value in iteritems(gene.__dict__):
                if attr not in do_not_copy_by_ref:
                    new_gene.__dict__[attr] = copy(
                        value) if attr == "formula" else value
            new_gene._model = new
            new.genes.append(new_gene)

        new.reactions = DictList()
        do_not_copy_by_ref = {"_model", "_metabolites", "_genes"}
        for reaction in self.reactions:
            new_reaction = reaction.__class__()
            for attr, value in iteritems(reaction.__dict__):
                if attr not in do_not_copy_by_ref:
                    new_reaction.__dict__[attr] = copy(value)
            new_reaction._model = new
            new.reactions.append(new_reaction)
            # update awareness
            for metabolite, stoic in iteritems(reaction._metabolites):
                new_met = new.metabolites.get_by_id(metabolite.id)
                new_reaction._metabolites[new_met] = stoic
                new_met._reaction.add(new_reaction)
            for gene in reaction._genes:
                new_gene = new.genes.get_by_id(gene.id)
                new_reaction._genes.add(new_gene)
                new_gene._reaction.add(new_reaction)

        for reaction in new.reactions:
            reaction._reset_var_cache()
        try:
            new._solver = deepcopy(self.solver)
            # Cplex has an issue with deep copies
        except Exception:  # pragma: no cover
            new._solver = copy(self.solver)  # pragma: no cover

        return new

    def add_metabolites(self, metabolite_list):
        """Will add a list of metabolites to the model object and add new
        constraints accordingly.

        The change is reverted upon exit when using the model as a context.

        Parameters
        ----------
        metabolite_list : A list of `cobra.core.Metabolite` objects

        """
        if not hasattr(metabolite_list, '__iter__'):
            metabolite_list = [metabolite_list]
        if len(metabolite_list) == 0:
            return None

        # First check whether the metabolites exist in the model
        metabolite_list = [x for x in metabolite_list
                           if x.id not in self.metabolites]
        for x in metabolite_list:
            x._model = self
        self.metabolites += metabolite_list

        # from cameo ...
        to_add = []
        for met in metabolite_list:
            if met.id not in self.constraints:
                constraint = self.problem.Constraint(
                    S.Zero, name=met.id, lb=0, ub=0)
                to_add += [constraint]

        self.add_cons_vars(to_add)

        context = get_context(self)
        if context:
            context(partial(self.metabolites.__isub__, metabolite_list))
            for x in metabolite_list:
                # Do we care?
                context(partial(setattr, x, '_model', None))

    def remove_metabolites(self, metabolite_list, destructive=False):
        """Remove a list of metabolites from the the object.

        The change is reverted upon exit when using the model as a context.

        Parameters
        ----------
        metabolite_list : list
            A list with `cobra.Metabolite` objects as elements.

        destructive : bool
            If False then the metabolite is removed from all
            associated reactions.  If True then all associated
            reactions are removed from the Model.

        """
        if not hasattr(metabolite_list, '__iter__'):
            metabolite_list = [metabolite_list]
        # Make sure metabolites exist in model
        metabolite_list = [x for x in metabolite_list
                           if x.id in self.metabolites]
        for x in metabolite_list:
            x._model = None

            if not destructive:
                for the_reaction in list(x._reaction):
                    the_coefficient = the_reaction._metabolites[x]
                    the_reaction.subtract_metabolites({x: the_coefficient})

            else:
                for x in list(x._reaction):
                    x.remove_from_model()

        self.metabolites -= metabolite_list

        to_remove = [self.solver.constraints[m.id] for m in metabolite_list]
        self.remove_cons_vars(to_remove)

        context = get_context(self)
        if context:
            context(partial(self.metabolites.__iadd__, metabolite_list))
            for x in metabolite_list:
                context(partial(setattr, x, '_model', self))

    def add_reaction(self, reaction):
        """Will add a cobra.Reaction object to the model, if
        reaction.id is not in self.reactions.

        Parameters
        ----------
        reaction : cobra.Reaction
            The reaction to add

        Deprecated (0.6). Use `~cobra.Model.add_reactions` instead
        """
        warn("add_reaction deprecated. Use add_reactions instead",
             DeprecationWarning)

        self.add_reactions([reaction])

    def add_reactions(self, reaction_list):
        """Will add a cobra.Reaction object to the model, if
        reaction.id is not in self.reactions.

        The change is reverted upon exit when using the model as a context.

        Parameters
        ----------
        reaction_list : list
            A list of `cobra.Reaction` objects

        """

        try:
            reaction_list = DictList(reaction_list)
        except TypeError:
            reaction_list = DictList([reaction_list])

        # Only add the reaction if one with the same ID is not already
        # present in the model.
        reactions_in_model = [i.id for i in reaction_list if i.id
                              in self.reactions]

        if len(reactions_in_model) > 0:
            raise Exception("Reactions already in the model: " +
                            ", ".join(reactions_in_model))

        context = get_context(self)

        # Add reactions. Also take care of genes and metabolites in the loop
        for reaction in reaction_list:
            reaction._reset_var_cache()
            reaction._model = self  # the reaction now points to the model
            # keys() is necessary because the dict will be modified during
            # the loop
            for metabolite in list(reaction._metabolites.keys()):
                # if the metabolite is not in the model, add it
                # should we be adding a copy instead.
                if metabolite not in self.metabolites:
                    self.add_metabolites(metabolite)
                # A copy of the metabolite exists in the model, the reaction
                # needs to point to the metabolite in the model.
                else:
                    stoichiometry = reaction._metabolites.pop(metabolite)
                    model_metabolite = self.metabolites.get_by_id(
                        metabolite.id)
                    reaction._metabolites[model_metabolite] = stoichiometry
                    model_metabolite._reaction.add(reaction)
                    if context:
                        context(partial(
                            model_metabolite._reaction.remove, reaction))

            for gene in list(reaction._genes):
                # If the gene is not in the model, add it
                if not self.genes.has_id(gene.id):
                    self.genes += [gene]
                    gene._model = self

                    if context:
                        # Remove the gene later
                        context(partial(self.genes.__isub__, [gene]))
                        context(partial(setattr, gene, '_model', None))

                # Otherwise, make the gene point to the one in the model
                else:
                    model_gene = self.genes.get_by_id(gene.id)
                    if model_gene is not gene:
                        reaction._dissociate_gene(gene)
                        reaction._associate_gene(model_gene)

        self.reactions += reaction_list

        if context:
            context(partial(self.reactions.__isub__, reaction_list))

        # from cameo ...
        self._populate_solver(reaction_list)

    def remove_reactions(self, reactions, delete=True,
                         remove_orphans=False):
        """Remove reactions from the model.

        The change is reverted upon exit when using the model as a context.

        Parameters
        ----------
        reactions : list
            A list with reactions (`cobra.Reaction`), or their id's, to remove

        delete : bool
            Whether or not the reactions should be deleted after removal.
            If the reactions are not deleted, those objects will be
            recreated with new metabolite and gene objects.

        remove_orphans : bool
            Remove orphaned genes and metabolites from the model as well

        """
        if isinstance(reactions, string_types) or hasattr(reactions, "id"):
            warn("need to pass in a list")
            reactions = [reactions]

        context = get_context(self)

        for reaction in reactions:
            try:
                reaction = self.reactions[self.reactions.index(reaction)]
            except ValueError:
                warn('%s not in %s' % (reaction, self))
            else:
                forward = reaction.forward_variable
                reverse = reaction.reverse_variable
                self.remove_cons_vars([forward, reverse])
                self.reactions.remove(reaction)
                reaction._model = None

                if context:
                    context(partial(setattr, reaction, '_model', self))
                    context(partial(self.reactions.add, reaction))

                for x in reaction._metabolites:
                    if reaction in x._reaction:
                        x._reaction.remove(reaction)
                        if context:
                            context(partial(x._reaction.add, reaction))
                        if remove_orphans and len(x._reaction) == 0:
                            self.remove_metabolites(x)

                for x in reaction._genes:
                    if reaction in x._reaction:
                        x._reaction.remove(reaction)
                        if context:
                            context(partial(x._reaction.add, reaction))

                        if remove_orphans and len(x._reaction) == 0:
                            self.genes.remove(x)
                            if context:
                                context(partial(self.genes.add, x))

                reaction._metabolites = {}
                reaction._genes = set()

    def add_cons_vars(self, what, **kwargs):
        """Add constraints and variables to the model's mathematical problem.

        Useful for variables and constraints that can not be expressed with
        reactions and simple lower and upper bounds.

        Additions are reversed upon exit if the model itself is used as
        context.

        Parameters
        ----------
        what : list or tuple of optlang variables or constraints.
           The variables or constraints to add to the model. Must be of
           class `optlang.interface.Variable` or
           `optlang.interface.Constraint`.
        **kwargs : keyword arguments
           Passed to solver.add()
        """
        add_cons_vars_to_problem(self, what, **kwargs)

    def remove_cons_vars(self, what):
        """Remove variables and constraints from the model's mathematical
        problem.

        Remove variables and constraints that were added directly to the
        model's underlying mathematical problem. Removals are reversed
        upon exit if the model itself is used as context.

        Parameters
        ----------
        what : list or tuple of optlang variables or constraints.
           The variables or constraints to add to the model. Must be of
           class `optlang.interface.Variable` or
           `optlang.interface.Constraint`.
        """
        remove_cons_vars_from_problem(self, what)

    @property
    def problem(self):
        """The interface to the model's underlying mathematical problem.

        Solutions to cobra models are obtained by formulating a mathematical
        problem and solving it. Cobrapy uses the optlang package to
        accomplish that and with this property you can get access to the
        problem interface directly.

        Returns
        -------
        optlang.interface
            The problem interface that defines methods for interacting with
            the problem and associated solver directly.
        """
        return self.solver.interface

    @property
    def variables(self):
        """The mathematical variables in the cobra model.

        In a cobra model, most variables are reactions. However,
        for specific use cases, it may also be useful to have other types of
        variables. This property defines all variables currently associated
        with the model's problem.

        Returns
        -------
        optlang.container.Container
            A container with all associated variables.
        """
        return self.solver.variables

    @property
    def constraints(self):
        """The constraints in the cobra model.

        In a cobra model, most constraints are metabolites and their
        stoichiometries. However, for specific use cases, it may also be
        useful to have other types of constraints. This property defines all
        constraints currently associated with the model's problem.

        Returns
        -------
        optlang.container.Container
            A container with all associated constraints.
        """
        return self.solver.constraints

    @property
    def exchanges(self):
        """Exchange reactions in model.

        Reactions that either don't have products or substrates.
        """
        return [rxn for rxn in self.reactions if rxn.boundary]

    def _populate_solver(self, reaction_list, metabolite_list=None):
        """Populate attached solver with constraints and variables that
        model the provided reactions.
        """
        constraint_terms = AutoVivification()
        to_add = []
        if metabolite_list is not None:
            for met in metabolite_list:
                to_add += [self.problem.Constraint(
                    S.Zero, name=met.id, lb=0, ub=0)]
        self.add_cons_vars(to_add)

        for reaction in reaction_list:

            reverse_lb, reverse_ub, forward_lb, forward_ub = \
                separate_forward_and_reverse_bounds(*reaction.bounds)

            forward_variable = self.problem.Variable(
                reaction.id, lb=forward_lb, ub=forward_ub)
            reverse_variable = self.problem.Variable(
                reaction.reverse_id, lb=reverse_lb, ub=reverse_ub)

            self.add_cons_vars([forward_variable, reverse_variable])
            self.solver.update()

            for metabolite, coeff in six.iteritems(reaction.metabolites):
                if metabolite.id in self.constraints:
                    constraint = self.constraints[metabolite.id]
                else:
                    constraint = self.problem.Constraint(
                        S.Zero,
                        name=metabolite.id,
                        lb=0, ub=0)
                    self.add_cons_vars(constraint, sloppy=True)

                constraint_terms[constraint][forward_variable] = coeff
                constraint_terms[constraint][reverse_variable] = -coeff

        self.solver.update()
        for constraint, terms in six.iteritems(constraint_terms):
            constraint.set_linear_coefficients(terms)

    def to_array_based_model(self, deepcopy_model=False, **kwargs):
        """Makes a `cobra.core.ArrayBasedModel` from a cobra.Model
        which may be used to perform linear algebra operations with the
        stoichiometric matrix.

        Deprecated (0.6). Use `cobra.util.array.create_stoichiometric_matrix`
        instead.

        Parameters
        ----------
        deepcopy_model : bool
            If False then the ArrayBasedModel points to the Model

        """
        warn("to_array_based_model is deprecated. "
             "use cobra.util.array.create_stoichiometric_matrix instead",
             DeprecationWarning)
        from cobra.core.arraybasedmodel import ArrayBasedModel
        return ArrayBasedModel(self, deepcopy_model=deepcopy_model, **kwargs)

    def optimize(self, objective_sense=None, **kwargs):
        """
        Optimize the model using flux balance analysis.

        Parameters
        ----------
        objective_sense : {None, 'maximize' 'minimize'}, optional
            Whether fluxes should be maximized or minimized. In case of None,
            the previous direction is used.
        solver : {None, 'glpk', 'cglpk', 'gurobi', 'cplex'}, optional
            If unspecified will use the currently defined `self.solver`
            otherwise it will use the given solver and update the attribute.
        quadratic_component : {None, scipy.sparse.dok_matrix}, optional
            The dimensions should be (n, n) where n is the number of
            reactions. This sets the quadratic component (Q) of the
            objective coefficient, adding :math:`\\frac{1}{2} v^T \cdot Q
            \cdot v` to the objective. Ignored for optlang based solvers.
        tolerance_feasibility : float
            Solver tolerance for feasibility. Ignored for optlang based
            solvers
        tolerance_markowitz : float
            Solver threshold during pivot. Ignored for optlang based solvers
        time_limit : float
            Maximum solver time (in seconds). Ignored for optlang based solvers

        Notes
        -----
        Only the most commonly used parameters are presented here.  Additional
        parameters for cobra.solvers may be available and specified with the
        appropriate keyword argument.

        """
        legacy, solver = choose_solver(self, solver=kwargs.get("solver"))
        original_direction = self.objective.direction

        if legacy:
            if objective_sense is None:
                objective_sense = {
                    "max": "maximize", "min": "minimize"}[original_direction]
            solution = optimize(self, objective_sense=objective_sense,
                                **kwargs)
            check_solver_status(solution.status)
            return solution

        self.solver = solver
        self.objective.direction = \
            {"maximize": "max", "minimize": "min"}.get(
                objective_sense, original_direction)
        self.solver.optimize()
        solution = get_solution(self)
        self.objective.direction = original_direction
        return solution

    def repair(self, rebuild_index=True, rebuild_relationships=True):
        """Update all indexes and pointers in a model

        Parameters
        ----------
        rebuild_index : bool
            rebuild the indices kept in reactions, metabolites and genes
        rebuild_relationships : bool
             reset all associations between genes, metabolites, model and
             then re-add them.
        """
        if rebuild_index:  # DictList indexes
            self.reactions._generate_index()
            self.metabolites._generate_index()
            self.genes._generate_index()
        if rebuild_relationships:
            for met in self.metabolites:
                met._reaction.clear()
            for gene in self.genes:
                gene._reaction.clear()
            for rxn in self.reactions:
                for met in rxn._metabolites:
                    met._reaction.add(rxn)
                for gene in rxn._genes:
                    gene._reaction.add(rxn)
        # point _model to self
        for l in (self.reactions, self.genes, self.metabolites):
            for e in l:
                e._model = self

    @property
    def objective(self):
        """Get or set the solver objective

        Before introduction of the optlang based problems,
        this function returned the objective reactions as a list. With
        optlang, the objective is not limited a simple linear summation of
        individual reaction fluxes, making that return value ambiguous.
        Henceforth, use `cobra.util.solver.linear_reaction_coefficients` to
        get a dictionary of reactions with their linear coefficients (empty
        if there are none)

        The set value can be dictionary (reactions as keys, linear
        coefficients as values), string (reaction identifier), int (reaction
        index), Reaction or problem.Objective or sympy expression
        directly interpreted as objectives.

        When using a `HistoryManager` context, this attribute can be set
        temporarily, reversed when the exiting the context.
        """
        return self.solver.objective

    @objective.setter
    def objective(self, value):
        if isinstance(value, sympy.Basic):
            value = self.problem.Objective(value, sloppy=False)
        if not isinstance(value, (dict, optlang.interface.Objective)):
            try:
                reactions = self.reactions.get_by_any(value)
            except KeyError:
                raise ValueError('invalid objective')
            value = {rxn: 1 for rxn in reactions}
        set_objective(self, value, additive=False)

    def summary(self, threshold=1E-8, fva=None, floatfmt='.3g', **kwargs):
        """Print a summary of the input and output fluxes of the model. This
        method requires the model to have been previously solved.

        Parameters
        ----------
        threshold : float
            tolerance for determining if a flux is zero (not printed)

        fva : int or None
            Whether or not to calculate and report flux variability in the
            output summary

        floatfmt : string
            format method for floats, passed to tabulate. Default is '.3g'.

        """
        from cobra.flux_analysis.summary import model_summary
        return model_summary(self, threshold=threshold, fva=fva,
                             floatfmt=floatfmt, **kwargs)

    def __enter__(self):
        """Record all future changes to the model, undoing them when a call to
        __exit__ is received"""

        # Create a new context and add it to the stack
        try:
            self._contexts.append(HistoryManager())
        except AttributeError:
            self._contexts = [HistoryManager()]

        return self

    def __exit__(self, type, value, traceback):
        """Pop the top context manager and trigger the undo functions"""
        context = self._contexts.pop()
        context.reset()
