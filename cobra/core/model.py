# -*- coding: utf-8 -*-

from __future__ import absolute_import

import logging
import types
from copy import copy, deepcopy
from functools import partial
from warnings import warn

import optlang
import six
from optlang.symbolics import Basic, Zero
from six import iteritems, string_types

from cobra.core.configuration import Configuration
from cobra.core.dictlist import DictList
from cobra.core.gene import Gene
from cobra.core.group import Group
from cobra.core.metabolite import Metabolite
from cobra.core.object import Object
from cobra.core.reaction import Reaction
from cobra.core.solution import get_solution
from cobra.exceptions import SolverNotFound
from cobra.medium import (
    find_boundary_types, find_external_compartment, sbo_terms)
from cobra.util.context import HistoryManager, get_context, resettable
from cobra.util.solver import (
    add_cons_vars_to_problem, assert_optimal, interface_to_str,
    remove_cons_vars_from_problem, set_objective, solvers)
from cobra.util.util import AutoVivification, format_long_string


LOGGER = logging.getLogger(__name__)
CONFIGURATION = Configuration()


class Model(Object):
    """Class representation for a cobra model

    Parameters
    ----------
    id_or_model : Model, string
        Either an existing Model object in which case a new model object is
        instantiated with the same properties as the original model,
        or an identifier to associate with the model as a string.
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
    groups : DictList
        A DictList where the key is the group identifier and the value a
        Group
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
            self.groups = DictList()  # A list of cobra.Groups
            # genes based on their ids {Gene.id: Gene}
            self._compartments = {}
            self._contexts = []

            # from cameo ...

            # if not hasattr(self, '_solver'):  # backwards compatibility
            # with older cobrapy pickles?

            interface = CONFIGURATION.solver
            self._solver = interface.Model()
            self._solver.objective = interface.Objective(Zero)
            self._populate_solver(self.reactions, self.metabolites)

            self._tolerance = None
            self.tolerance = CONFIGURATION.tolerance

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
            '%s is not a valid solver interface. Pick from %s.' % (
                value, list(solvers)))
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
        self._solver = interface.Model.clone(self._solver)

    @property
    def tolerance(self):
        return self._tolerance

    @tolerance.setter
    def tolerance(self, value):
        solver_tolerances = self._solver.configuration.tolerances

        try:
            solver_tolerances.feasibility = value
        except AttributeError:
            LOGGER.info("The current solver doesn't allow setting"
                        "feasibility tolerance.")

        try:
            solver_tolerances.optimality = value
        except AttributeError:
            LOGGER.info("The current solver doesn't allow setting"
                        "optimality tolerance.")

        try:
            solver_tolerances.integrality = value
        except AttributeError:
            LOGGER.info("The current solver doesn't allow setting"
                        "integrality tolerance.")

        self._tolerance = value

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
        warn('use Model.compartments instead', DeprecationWarning)
        return {met.compartment for met in self.metabolites
                if met.compartment is not None}

    @property
    def compartments(self):
        return {met.compartment: self._compartments.get(met.compartment, '')
                for met in self.metabolites if met.compartment is not None}

    @compartments.setter
    def compartments(self, value):
        """Get or set the dictionary of current compartment descriptions.

        Assigning a dictionary to this property updates the model's
        dictionary of compartment descriptions with the new values.

        Parameters
        ----------
        value : dict
            Dictionary mapping compartments abbreviations to full names.

        Examples
        --------
        >>> import cobra.test
        >>> model = cobra.test.create_test_model("textbook")
        >>> model.compartments = {'c': 'the cytosol'}
        {'c': 'the cytosol', 'e': 'extracellular'}
        """
        self._compartments.update(value)

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
        exchange_rxns = frozenset(self.exchanges)
        for rxn_id, bound in iteritems(medium):
            rxn = self.reactions.get_by_id(rxn_id)
            if rxn not in exchange_rxns:
                LOGGER.warn("%s does not seem to be an"
                            " an exchange reaction. Applying bounds anyway.",
                            rxn.id)
            media_rxns.append(rxn)
            set_active_bound(rxn, bound)

        media_rxns = frozenset(media_rxns)

        # Turn off reactions not present in media
        for rxn in (exchange_rxns - media_rxns):
            set_active_bound(rxn, 0)

    def __add__(self, other_model):
        """Add the content of another model to this model (+).

        The model is copied as a new object, with a new model identifier,
        and copies of all the reactions in the other model are added to this
        model. The objective is the sum of the objective expressions for the
        two models.
        """
        warn('use model.merge instead', DeprecationWarning)
        return self.merge(other_model, objective='sum', inplace=False)

    def __iadd__(self, other_model):
        """Incrementally add the content of another model to this model (+=).

        Copies of all the reactions in the other model are added to this
        model. The objective is the sum of the objective expressions for the
        two models.
        """
        warn('use model.merge instead', DeprecationWarning)
        return self.merge(other_model, objective='sum', inplace=True)

    def copy(self):
        """Provides a partial 'deepcopy' of the Model.  All of the Metabolite,
        Gene, and Reaction objects are created anew but in a faster fashion
        than deepcopy
        """
        new = self.__class__()
        do_not_copy_by_ref = {"metabolites", "reactions", "genes", "notes",
                              "annotation", "groups"}
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

        new.groups = DictList()
        do_not_copy_by_ref = {"_model", "_members"}
        # Groups can be members of other groups. We initialize them first and
        # then update their members.
        for group in self.groups:
            new_group = group.__class__(group.id)
            for attr, value in iteritems(group.__dict__):
                if attr not in do_not_copy_by_ref:
                    new_group.__dict__[attr] = copy(value)
            new_group._model = new
            new.groups.append(new_group)
        for group in self.groups:
            new_group = new.groups.get_by_id(group.id)
            # update awareness, as in the reaction copies
            new_objects = []
            for member in group.members:
                if isinstance(member, Metabolite):
                    new_object = new.metabolites.get_by_id(member.id)
                elif isinstance(member, Reaction):
                    new_object = new.reactions.get_by_id(member.id)
                elif isinstance(member, Gene):
                    new_object = new.genes.get_by_id(member.id)
                elif isinstance(member, Group):
                    new_object = new.genes.get_by_id(member.id)
                else:
                    raise TypeError(
                        "The group member {!r} is unexpectedly not a "
                        "metabolite, reaction, gene, nor another "
                        "group.".format(member))
                new_objects.append(new_object)
            new_group.add_members(new_objects)

        try:
            new._solver = deepcopy(self.solver)
            # Cplex has an issue with deep copies
        except Exception:  # pragma: no cover
            new._solver = copy(self.solver)  # pragma: no cover

        # it doesn't make sense to retain the context of a copied model so
        # assign a new empty context
        new._contexts = list()

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

        bad_ids = [m for m in metabolite_list
                   if not isinstance(m.id, string_types) or len(m.id) < 1]
        if len(bad_ids) != 0:
            raise ValueError('invalid identifiers in {}'.format(repr(bad_ids)))

        for x in metabolite_list:
            x._model = self
        self.metabolites += metabolite_list

        # from cameo ...
        to_add = []
        for met in metabolite_list:
            if met.id not in self.constraints:
                constraint = self.problem.Constraint(
                    Zero, name=met.id, lb=0, ub=0)
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

            # remove reference to the metabolite in all groups
            associated_groups = self.get_associated_groups(x)
            for group in associated_groups:
                group.remove_members(x)

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

    def add_boundary(self, metabolite, type="exchange", reaction_id=None,
                     lb=None, ub=None, sbo_term=None):
        """
        Add a boundary reaction for a given metabolite.

        There are three different types of pre-defined boundary reactions:
        exchange, demand, and sink reactions.
        An exchange reaction is a reversible, unbalanced reaction that adds
        to or removes an extracellular metabolite from the extracellular
        compartment.
        A demand reaction is an irreversible reaction that consumes an
        intracellular metabolite.
        A sink is similar to an exchange but specifically for intracellular
        metabolites.

        If you set the reaction `type` to something else, you must specify the
        desired identifier of the created reaction along with its upper and
        lower bound. The name will be given by the metabolite name and the
        given `type`.

        Parameters
        ----------
        metabolite : cobra.Metabolite
            Any given metabolite. The compartment is not checked but you are
            encouraged to stick to the definition of exchanges and sinks.
        type : str, {"exchange", "demand", "sink"}
            Using one of the pre-defined reaction types is easiest. If you
            want to create your own kind of boundary reaction choose
            any other string, e.g., 'my-boundary'.
        reaction_id : str, optional
            The ID of the resulting reaction. This takes precedence over the
            auto-generated identifiers but beware that it might make boundary
            reactions harder to identify afterwards when using `model.boundary`
            or specifically `model.exchanges` etc.
        lb : float, optional
            The lower bound of the resulting reaction.
        ub : float, optional
            The upper bound of the resulting reaction.
        sbo_term : str, optional
            A correct SBO term is set for the available types. If a custom
            type is chosen, a suitable SBO term should also be set.

        Returns
        -------
        cobra.Reaction
            The created boundary reaction.

        Examples
        --------
        >>> import cobra.test
        >>> model = cobra.test.create_test_model("textbook")
        >>> demand = model.add_boundary(model.metabolites.atp_c, type="demand")
        >>> demand.id
        'DM_atp_c'
        >>> demand.name
        'ATP demand'
        >>> demand.bounds
        (0, 1000.0)
        >>> demand.build_reaction_string()
        'atp_c --> '

        """
        ub = CONFIGURATION.upper_bound if ub is None else ub
        lb = CONFIGURATION.lower_bound if lb is None else lb
        types = {
            "exchange": ("EX", lb, ub, sbo_terms["exchange"]),
            "demand": ("DM", 0, ub, sbo_terms["demand"]),
            "sink": ("SK", lb, ub, sbo_terms["sink"])
        }
        if type == "exchange":
            external = find_external_compartment(self)
            if metabolite.compartment != external:
                raise ValueError("The metabolite is not an external metabolite"
                                 " (compartment is `%s` but should be `%s`). "
                                 "Did you mean to add a demand or sink? "
                                 "If not, either change its compartment or "
                                 "rename the model compartments to fix this." %
                                 (metabolite.compartment, external))
        if type in types:
            prefix, lb, ub, default_term = types[type]
            if reaction_id is None:
                reaction_id = "{}_{}".format(prefix, metabolite.id)
            if sbo_term is None:
                sbo_term = default_term
        if reaction_id is None:
            raise ValueError(
                "Custom types of boundary reactions require a custom "
                "identifier. Please set the `reaction_id`.")
        if reaction_id in self.reactions:
            raise ValueError(
                "Boundary reaction '{}' already exists.".format(reaction_id))
        name = "{} {}".format(metabolite.name, type)
        rxn = Reaction(id=reaction_id, name=name, lower_bound=lb,
                       upper_bound=ub)
        rxn.add_metabolites({metabolite: -1})
        if sbo_term:
            rxn.annotation["sbo"] = sbo_term
        self.add_reactions([rxn])
        return rxn

    def add_reactions(self, reaction_list):
        """Add reactions to the model.

        Reactions with identifiers identical to a reaction already in the
        model are ignored.

        The change is reverted upon exit when using the model as a context.

        Parameters
        ----------
        reaction_list : list
            A list of `cobra.Reaction` objects
        """
        def existing_filter(rxn):
            if rxn.id in self.reactions:
                LOGGER.warning(
                    "Ignoring reaction '%s' since it already exists.", rxn.id)
                return False
            return True

        # First check whether the reactions exist in the model.
        pruned = DictList(filter(existing_filter, reaction_list))

        context = get_context(self)

        # Add reactions. Also take care of genes and metabolites in the loop.
        for reaction in pruned:
            reaction._model = self
            # Build a `list()` because the dict will be modified in the loop.
            for metabolite in list(reaction.metabolites):
                # TODO: Should we add a copy of the metabolite instead?
                if metabolite not in self.metabolites:
                    self.add_metabolites(metabolite)
                # A copy of the metabolite exists in the model, the reaction
                # needs to point to the metabolite in the model.
                else:
                    # FIXME: Modifying 'private' attributes is horrible.
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

        self.reactions += pruned

        if context:
            context(partial(self.reactions.__isub__, pruned))

        # from cameo ...
        self._populate_solver(pruned)

    def remove_reactions(self, reactions, remove_orphans=False):
        """Remove reactions from the model.

        The change is reverted upon exit when using the model as a context.

        Parameters
        ----------
        reactions : list
            A list with reactions (`cobra.Reaction`), or their id's, to remove

        remove_orphans : bool
            Remove orphaned genes and metabolites from the model as well

        """
        if isinstance(reactions, string_types) or hasattr(reactions, "id"):
            warn("need to pass in a list")
            reactions = [reactions]

        context = get_context(self)

        for reaction in reactions:

            # Make sure the reaction is in the model
            try:
                reaction = self.reactions[self.reactions.index(reaction)]
            except ValueError:
                warn('%s not in %s' % (reaction, self))

            else:
                forward = reaction.forward_variable
                reverse = reaction.reverse_variable

                if context:

                    obj_coef = reaction.objective_coefficient

                    if obj_coef != 0:
                        context(partial(
                            self.solver.objective.set_linear_coefficients,
                            {forward: obj_coef, reverse: -obj_coef}))

                    context(partial(self._populate_solver, [reaction]))
                    context(partial(setattr, reaction, '_model', self))
                    context(partial(self.reactions.add, reaction))

                self.remove_cons_vars([forward, reverse])
                self.reactions.remove(reaction)
                reaction._model = None

                for met in reaction._metabolites:
                    if reaction in met._reaction:
                        met._reaction.remove(reaction)
                        if context:
                            context(partial(met._reaction.add, reaction))
                        if remove_orphans and len(met._reaction) == 0:
                            self.remove_metabolites(met)

                for gene in reaction._genes:
                    if reaction in gene._reaction:
                        gene._reaction.remove(reaction)
                        if context:
                            context(partial(gene._reaction.add, reaction))

                        if remove_orphans and len(gene._reaction) == 0:
                            self.genes.remove(gene)
                            if context:
                                context(partial(self.genes.add, gene))

                # remove reference to the reaction in all groups
                associated_groups = self.get_associated_groups(reaction)
                for group in associated_groups:
                    group.remove_members(reaction)

    def add_groups(self, group_list):
        """Add groups to the model.

        Groups with identifiers identical to a group already in the model are
        ignored.

        If any group contains members that are not in the model, these members
        are added to the model as well. Only metabolites, reactions, and genes
        can have groups.

        Parameters
        ----------
        group_list : list
            A list of `cobra.Group` objects to add to the model.
        """

        def existing_filter(group):
            if group.id in self.groups:
                LOGGER.warning(
                    "Ignoring group '%s' since it already exists.", group.id)
                return False
            return True

        if isinstance(group_list, string_types) or \
                hasattr(group_list, "id"):
            warn("need to pass in a list")
            group_list = [group_list]

        pruned = DictList(filter(existing_filter, group_list))

        for group in pruned:
            group._model = self
            for member in group.members:
                # If the member is not associated with the model, add it
                if isinstance(member, Metabolite):
                    if member not in self.metabolites:
                        self.add_metabolites([member])
                if isinstance(member, Reaction):
                    if member not in self.reactions:
                        self.add_reactions([member])
                # TODO(midnighter): `add_genes` method does not exist.
                # if isinstance(member, Gene):
                #     if member not in self.genes:
                #         self.add_genes([member])

            self.groups += [group]

    def remove_groups(self, group_list):
        """Remove groups from the model.

        Members of each group are not removed
        from the model (i.e. metabolites, reactions, and genes in the group
        stay in the model after any groups containing them are removed).

        Parameters
        ----------
        group_list : list
            A list of `cobra.Group` objects to remove from the model.
        """

        if isinstance(group_list, string_types) or \
                hasattr(group_list, "id"):
            warn("need to pass in a list")
            group_list = [group_list]

        for group in group_list:
            # make sure the group is in the model
            if group.id not in self.groups:
                LOGGER.warning("%r not in %r. Ignored.", group, self)
            else:
                self.groups.remove(group)
                group._model = None

    def get_associated_groups(self, element):
        """Returns a list of groups that an element (reaction, metabolite, gene)
        is associated with.

        Parameters
        ----------
        element: `cobra.Reaction`, `cobra.Metabolite`, or `cobra.Gene`

        Returns
        -------
        list of `cobra.Group`
            All groups that the provided object is a member of
        """
        # check whether the element is associated with the model
        return [g for g in self.groups if element in g.members]

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
    def boundary(self):
        """Boundary reactions in the model.
        Reactions that either have no substrate or product.
        """
        return [rxn for rxn in self.reactions if rxn.boundary]

    @property
    def exchanges(self):
        """Exchange reactions in model.
        Reactions that exchange mass with the exterior. Uses annotations
        and heuristics to exclude non-exchanges such as sink reactions.
        """
        return find_boundary_types(self, "exchange", None)

    @property
    def demands(self):
        """Demand reactions in model.
        Irreversible reactions that accumulate or consume a metabolite in
        the inside of the model.
        """
        return find_boundary_types(self, "demand", None)

    @property
    def sinks(self):
        """Sink reactions in model.
        Reversible reactions that accumulate or consume a metabolite in
        the inside of the model.
        """
        return find_boundary_types(self, "sink", None)

    def _populate_solver(self, reaction_list, metabolite_list=None):
        """Populate attached solver with constraints and variables that
        model the provided reactions.
        """
        constraint_terms = AutoVivification()
        to_add = []
        if metabolite_list is not None:
            for met in metabolite_list:
                to_add += [self.problem.Constraint(
                    Zero, name=met.id, lb=0, ub=0)]
        self.add_cons_vars(to_add)

        for reaction in reaction_list:
            if reaction.id not in self.variables:
                forward_variable = self.problem.Variable(reaction.id)
                reverse_variable = self.problem.Variable(reaction.reverse_id)
                self.add_cons_vars([forward_variable, reverse_variable])
            else:
                reaction = self.reactions.get_by_id(reaction.id)
                forward_variable = reaction.forward_variable
                reverse_variable = reaction.reverse_variable
            for metabolite, coeff in six.iteritems(reaction.metabolites):
                if metabolite.id in self.constraints:
                    constraint = self.constraints[metabolite.id]
                else:
                    constraint = self.problem.Constraint(
                        Zero,
                        name=metabolite.id,
                        lb=0, ub=0)
                    self.add_cons_vars(constraint, sloppy=True)
                constraint_terms[constraint][forward_variable] = coeff
                constraint_terms[constraint][reverse_variable] = -coeff

        self.solver.update()
        for reaction in reaction_list:
            reaction = self.reactions.get_by_id(reaction.id)
            reaction.update_variable_bounds()
        for constraint, terms in six.iteritems(constraint_terms):
            constraint.set_linear_coefficients(terms)

    def slim_optimize(self, error_value=float('nan'), message=None):
        """Optimize model without creating a solution object.

        Creating a full solution object implies fetching shadow prices and
        flux values for all reactions and metabolites from the solver
        object. This necessarily takes some time and in cases where only one
        or two values are of interest, it is recommended to instead use this
        function which does not create a solution object returning only the
        value of the objective. Note however that the `optimize()` function
        uses efficient means to fetch values so if you need fluxes/shadow
        prices for more than say 4 reactions/metabolites, then the total
        speed increase of `slim_optimize` versus `optimize` is  expected to
        be small or even negative depending on how you fetch the values
        after optimization.

        Parameters
        ----------
        error_value : float, None
           The value to return if optimization failed due to e.g.
           infeasibility. If None, raise `OptimizationError` if the
           optimization fails.
        message : string
           Error message to use if the model optimization did not succeed.

        Returns
        -------
        float
            The objective value.
        """
        self.solver.optimize()
        if self.solver.status == optlang.interface.OPTIMAL:
            return self.solver.objective.value
        elif error_value is not None:
            return error_value
        else:
            assert_optimal(self, message)

    def optimize(self, objective_sense=None, raise_error=False):
        """
        Optimize the model using flux balance analysis.

        Parameters
        ----------
        objective_sense : {None, 'maximize' 'minimize'}, optional
            Whether fluxes should be maximized or minimized. In case of None,
            the previous direction is used.
        raise_error : bool
            If true, raise an OptimizationError if solver status is not
             optimal.

        Notes
        -----
        Only the most commonly used parameters are presented here.  Additional
        parameters for cobra.solvers may be available and specified with the
        appropriate keyword argument.

        """
        original_direction = self.objective.direction
        self.objective.direction = \
            {"maximize": "max", "minimize": "min"}.get(
                objective_sense, original_direction)
        self.slim_optimize()
        solution = get_solution(self, raise_error=raise_error)
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
            self.groups._generate_index()
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
        for l in (self.reactions, self.genes, self.metabolites, self.groups):
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
        if isinstance(value, Basic):
            value = self.problem.Objective(value, sloppy=False)
        if not isinstance(value, (dict, optlang.interface.Objective)):
            try:
                reactions = self.reactions.get_by_any(value)
            except KeyError:
                raise ValueError('invalid objective')
            value = {rxn: 1 for rxn in reactions}
        set_objective(self, value, additive=False)

    @property
    def objective_direction(self):
        """
        Get or set the objective direction.

        When using a `HistoryManager` context, this attribute can be set
        temporarily, reversed when exiting the context.

        """
        return self.solver.objective.direction

    @objective_direction.setter
    @resettable
    def objective_direction(self, value):
        value = value.lower()
        if value.startswith("max"):
            self.solver.objective.direction = "max"
        elif value.startswith("min"):
            self.solver.objective.direction = "min"
        else:
            raise ValueError("Unknown objective direction '{}'.".format(value))

    def summary(self, solution=None, threshold=1E-06, fva=None, names=False,
                floatfmt='.3g'):
        """
        Print a summary of the input and output fluxes of the model.

        Parameters
        ----------
        solution: cobra.Solution, optional
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
        from cobra.flux_analysis.summary import model_summary
        return model_summary(self, solution=solution, threshold=threshold,
                             fva=fva, names=names, floatfmt=floatfmt)

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

    def merge(self, right, prefix_existing=None, inplace=True,
              objective='left'):
        """Merge two models to create a model with the reactions from both
        models.

        Custom constraints and variables from right models are also copied
        to left model, however note that, constraints and variables are
        assumed to be the same if they have the same name.

        right : cobra.Model
            The model to add reactions from
        prefix_existing : string
            Prefix the reaction identifier in the right that already exist
            in the left model with this string.
        inplace : bool
            Add reactions from right directly to left model object.
            Otherwise, create a new model leaving the left model untouched.
            When done within the model as context, changes to the models are
            reverted upon exit.
        objective : string
            One of 'left', 'right' or 'sum' for setting the objective of the
            resulting model to that of the corresponding model or the sum of
            both.
        """
        if inplace:
            new_model = self
        else:
            new_model = self.copy()
            new_model.id = '{}_{}'.format(self.id, right.id)
        new_reactions = deepcopy(right.reactions)
        if prefix_existing is not None:
            existing = new_reactions.query(
                lambda rxn: rxn.id in self.reactions)
            for reaction in existing:
                reaction.id = '{}{}'.format(prefix_existing, reaction.id)
        new_model.add_reactions(new_reactions)
        interface = new_model.problem
        new_vars = [interface.Variable.clone(v) for v in right.variables if
                    v.name not in new_model.variables]
        new_model.add_cons_vars(new_vars)
        new_cons = [interface.Constraint.clone(c, model=new_model.solver)
                    for c in right.constraints if
                    c.name not in new_model.constraints]
        new_model.add_cons_vars(new_cons, sloppy=True)
        new_model.objective = dict(
            left=self.objective,
            right=right.objective,
            sum=self.objective.expression + right.objective.expression
        )[objective]
        return new_model

    def _repr_html_(self):
        return """
        <table>
            <tr>
                <td><strong>Name</strong></td>
                <td>{name}</td>
            </tr><tr>
                <td><strong>Memory address</strong></td>
                <td>{address}</td>
            </tr><tr>
                <td><strong>Number of metabolites</strong></td>
                <td>{num_metabolites}</td>
            </tr><tr>
                <td><strong>Number of reactions</strong></td>
                <td>{num_reactions}</td>
            </tr><tr>
                <td><strong>Number of groups</strong></td>
                <td>{num_groups}</td>
            </tr><tr>
                <td><strong>Objective expression</strong></td>
                <td>{objective}</td>
            </tr><tr>
                <td><strong>Compartments</strong></td>
                <td>{compartments}</td>
            </tr>
          </table>""".format(
            name=self.id,
            address='0x0%x' % id(self),
            num_metabolites=len(self.metabolites),
            num_reactions=len(self.reactions),
            num_groups=len(self.groups),
            objective=format_long_string(str(self.objective.expression), 100),
            compartments=", ".join(
                v if v else k for k, v in iteritems(self.compartments)
            ))
