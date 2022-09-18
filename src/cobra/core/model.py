"""Define the Model class."""

import logging
from copy import copy, deepcopy
from functools import partial
from types import ModuleType
from typing import TYPE_CHECKING, Dict, Iterable, List, Optional, Tuple, Union
from warnings import warn

import optlang
from optlang.symbolics import Basic, Zero

from ..medium import find_boundary_types, find_external_compartment, sbo_terms
from ..util.context import HistoryManager, get_context, resettable
from ..util.solver import (
    add_cons_vars_to_problem,
    assert_optimal,
    check_solver,
    interface_to_str,
    remove_cons_vars_from_problem,
    set_objective,
)
from ..util.util import AutoVivification, format_long_string
from .configuration import Configuration
from .dictlist import DictList
from .gene import Gene
from .group import Group
from .metabolite import Metabolite
from .object import Object
from .reaction import Reaction
from .solution import get_solution


if TYPE_CHECKING:

    import pandas as pd
    from optlang.container import Container

    from cobra import Solution
    from cobra.summary import ModelSummary
    from cobra.util.solver import CONS_VARS

logger = logging.getLogger(__name__)
configuration = Configuration()


class Model(Object):
    """Class representation for a cobra model.

    Parameters
    ----------
    id_or_model: str or Model, optional
        String to use as model id, or actual model to base new model one.
        If string, it is used as id. If model, a new model object is
        instantiated with the same properties as the original model (default None).
    name: str, optional
        Human readable string to be model description (default None).

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
    """

    def __init__(
        self, id_or_model: Union[str, "Model", None] = None, name: Optional[str] = None
    ) -> None:
        """Initialize the Model."""
        if isinstance(id_or_model, Model):
            Object.__init__(self, name=name)
            self.__setstate__(id_or_model.__dict__)
            self._solver = id_or_model.solver
        else:
            Object.__init__(self, id_or_model, name=name)
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

            interface = check_solver(configuration.solver)
            self._solver = interface.Model()
            self._solver.objective = interface.Objective(Zero)
            self._populate_solver(self.reactions, self.metabolites)

            self._tolerance = None
            self.tolerance = configuration.tolerance

    def __setstate__(self, state: Dict) -> None:
        """Make sure all cobra.Objects in the model point to the model.

        Parameters
        ----------
        state: dict
        """
        self.__dict__.update(state)
        for y in ["reactions", "genes", "metabolites"]:
            for x in getattr(self, y):
                x._model = self
        if not hasattr(self, "name"):
            self.name = None

    def __getstate__(self) -> Dict:
        """Get state for serialization.

        Ensures that the context stack is cleared prior to serialization,
        since partial functions cannot be pickled reliably.

        Returns
        -------
        odict: Dict
            A dictionary of state, based on self.__dict__.
        """
        odict = self.__dict__.copy()
        odict["_contexts"] = []
        return odict

    @property
    def solver(self) -> "optlang.interface.Model":
        """Get the attached solver instance.

        The associated the solver object, which manages the interaction with
        the associated solver, e.g. glpk.

        This property is useful for accessing the optimization problem
        directly and to define additional non-metabolic constraints.

        Examples
        --------
        >>> from cobra.io import load_model
        >>> model = load_model("textbook")
        >>> new = model.problem.Constraint(model.objective.expression, lb=0.99)
        >>> model.solver.add(new)
        """
        return self._solver

    @solver.setter
    @resettable
    def solver(self, value: Union[str, ModuleType]) -> None:
        """Set the attached solver instance.

        The associated the solver object, which manages the interaction with
        the associated solver, e.g. glpk.

        This property is useful for accessing the optimization problem
        directly and to define additional non-metabolic constraints.

        Parameters
        ----------
        value: ModuleType or str
            The solver to set, which will be checked by `check_solver()`.
        """
        interface = check_solver(value)

        # Do nothing if the solver did not change
        if self.problem == interface:
            return
        self._solver = interface.Model.clone(self._solver)

    @property
    def tolerance(self) -> float:
        """Get the tolerance.

        Returns
        -------
        float
            The tolerance of the mdoel.
        """
        return self._tolerance

    @tolerance.setter
    def tolerance(self, value: float) -> None:
        """Set the tolerance.

        Parameters
        ----------
        value: float
            Value to set tolerance.
        """
        solver_tolerances = self.solver.configuration.tolerances

        try:
            solver_tolerances.feasibility = value
        except AttributeError:
            logger.info(
                f"The current solver interface {interface_to_str(self.problem)} "
                f"doesn't support setting the feasibility tolerance."
            )

        try:
            solver_tolerances.optimality = value
        except AttributeError:
            logger.info(
                f"The current solver interface {interface_to_str(self.problem)} "
                f"doesn't support setting the optimality tolerance."
            )

        try:
            solver_tolerances.integrality = value
        except AttributeError:
            logger.info(
                f"The current solver interface {interface_to_str(self.problem)} "
                f"doesn't support setting the integrality tolerance."
            )

        self._tolerance = value

    @property
    def compartments(self) -> Dict:
        """Return all metabolites' compartments.

        Returns
        -------
        dict
            A dictionary of metabolite compartments, where the keys are the short
            version (one letter version) of the compartmetns, and the values are the
            full names (if they exist).
        """
        return {
            met.compartment: self._compartments.get(met.compartment, "")
            for met in self.metabolites
            if met.compartment is not None
        }

    @compartments.setter
    def compartments(self, value: Dict) -> None:
        """Get or set the dictionary of current compartment descriptions.

        Assigning a dictionary to this property updates the model's
        dictionary of compartment descriptions with the new values.

        Parameters
        ----------
        value : dict
            Dictionary mapping compartments abbreviations to full names.

        Examples
        --------
        >>> from cobra.io import load_model
        >>> model = load_model("textbook")
        >>> model.compartments = {'c': 'the cytosol'}
        >>> model.compartments
        {'c': 'the cytosol', 'e': 'extracellular'}
        """
        self._compartments.update(value)

    @property
    def medium(self) -> Dict[str, float]:
        """Get the constraints on the model exchanges.

        `model.medium` returns a dictionary of the bounds for each of the
        boundary reactions, in the form of `{rxn_id: bound}`, where `bound`
        specifies the absolute value of the bound in direction of metabolite
        creation (i.e., lower_bound for `met <--`, upper_bound for `met -->`)

        Returns
        -------
        Dict[str, float]
            A dictionary with rxn.id (str) as key, bound (float) as value.
        """

        def is_active(reaction: Reaction) -> bool:
            """Determine if boundary reaction permits flux towards creating metabolites.

            Parameters
            ----------
            reaction: cobra.Reaction

            Returns
            -------
            bool
                True if reaction produces metaoblites and has upper_bound above 0
                or if reaction consumes metabolites and has lower_bound below 0 (so
                could be reversed).
            """
            return (bool(reaction.products) and (reaction.upper_bound > 0)) or (
                bool(reaction.reactants) and (reaction.lower_bound < 0)
            )

        def get_active_bound(reaction: Reaction) -> float:
            """For an active boundary reaction, return the relevant bound.

            Parameters
            ----------
            reaction: cobra.Reaction

            Returns
            -------
            float:
                upper or minus lower bound, depenending if the reaction produces or
                consumes metaoblties.
            """
            if reaction.reactants:
                return -reaction.lower_bound
            elif reaction.products:
                return reaction.upper_bound

        return {
            rxn.id: get_active_bound(rxn) for rxn in self.exchanges if is_active(rxn)
        }

    @medium.setter
    def medium(self, medium: Dict[str, float]) -> None:
        """Set the constraints on the model exchanges.

        `model.medium` returns a dictionary of the bounds for each of the
        boundary reactions, in the form of `{rxn_id: rxn_bound}`, where `rxn_bound`
        specifies the absolute value of the bound in direction of metabolite
        creation (i.e., lower_bound for `met <--`, upper_bound for `met -->`)

        Parameters
        ----------
        medium: dict
            The medium to initialize. medium should be a dictionary defining
            `{rxn_id: bound}` pairs.
        """

        def set_active_bound(reaction: Reaction, bound: float) -> None:
            """Set active bound.

            Parameters
            ----------
            reaction: cobra.Reaction
                Reaction to set
            bound: float
                Value to set bound to. The bound is reversed and set as lower bound
                if reaction has reactants (metabolites that are consumed). If reaction
                has reactants, it seems the upper bound won't be set.
            """
            if reaction.reactants:
                reaction.lower_bound = -bound
            elif reaction.products:
                reaction.upper_bound = bound

        # Set the given media bounds
        media_rxns = []
        exchange_rxns = frozenset(self.exchanges)
        for rxn_id, rxn_bound in medium.items():
            rxn = self.reactions.get_by_id(rxn_id)
            if rxn not in exchange_rxns:
                logger.warning(
                    f"{rxn.id} does not seem to be an an exchange reaction. "
                    f"Applying bounds anyway."
                )
            media_rxns.append(rxn)
            # noinspection PyTypeChecker
            set_active_bound(rxn, rxn_bound)

        frozen_media_rxns = frozenset(media_rxns)

        # Turn off reactions not present in media
        for rxn in exchange_rxns - frozen_media_rxns:
            is_export = rxn.reactants and not rxn.products
            set_active_bound(
                rxn, min(0.0, -rxn.lower_bound if is_export else rxn.upper_bound)
            )

    def copy(self) -> "Model":
        """Provide a partial 'deepcopy' of the Model.

        All the Metabolite, Gene, and Reaction objects are created anew but
        in a faster fashion than deepcopy.

        Returns
        -------
        cobra.Model: new model copy
        """
        new = self.__class__()
        do_not_copy_by_ref = {
            "metabolites",
            "reactions",
            "genes",
            "notes",
            "annotation",
            "groups",
        }
        for attr in self.__dict__:
            if attr not in do_not_copy_by_ref:
                new.__dict__[attr] = self.__dict__[attr]
        new.notes = deepcopy(self.notes)
        new.annotation = deepcopy(self.annotation)

        new.metabolites = DictList()
        do_not_copy_by_ref = {"_reaction", "_model"}
        for metabolite in self.metabolites:
            new_met = metabolite.__class__()
            for attr, value in metabolite.__dict__.items():
                if attr not in do_not_copy_by_ref:
                    new_met.__dict__[attr] = copy(value) if attr == "formula" else value
            new_met._model = new
            new.metabolites.append(new_met)

        new.genes = DictList()
        for gene in self.genes:
            new_gene = gene.__class__(None)
            for attr, value in gene.__dict__.items():
                if attr not in do_not_copy_by_ref:
                    new_gene.__dict__[attr] = (
                        copy(value) if attr == "formula" else value
                    )
            new_gene._model = new
            new.genes.append(new_gene)

        new.reactions = DictList()
        do_not_copy_by_ref = {"_model", "_metabolites", "_genes"}
        for reaction in self.reactions:
            new_reaction = reaction.__class__()
            for attr, value in reaction.__dict__.items():
                if attr not in do_not_copy_by_ref:
                    new_reaction.__dict__[attr] = copy(value)
            new_reaction._model = new
            new.reactions.append(new_reaction)
            # update awareness
            for metabolite, stoic in reaction._metabolites.items():
                new_met = new.metabolites.get_by_id(metabolite.id)
                new_reaction._metabolites[new_met] = stoic
                new_met._reaction.add(new_reaction)
            new_reaction.update_genes_from_gpr()

        new.groups = DictList()
        do_not_copy_by_ref = {"_model", "_members"}
        # Groups can be members of other groups. We initialize them first and
        # then update their members.
        for group in self.groups:
            new_group: Group = group.__class__(group.id)
            for attr, value in group.__dict__.items():
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
                    new_object = new.groups.get_by_id(member.id)
                else:
                    raise TypeError(
                        f"The group member {member!r} is unexpectedly not a "
                        f"metabolite, reaction, gene, nor another group."
                    )
                new_objects.append(new_object)
            new_group.add_members(new_objects)

        try:
            new._solver = deepcopy(self.solver)
            # Cplex has an issue with deep copies
        except Exception:  # pragma: no cover
            new._solver = copy(self.solver)  # pragma: no cover

        # it doesn't make sense to retain the context of a copied model so
        # assign a new empty context
        new._contexts = []

        return new

    def add_metabolites(self, metabolite_list: Union[List, Metabolite]) -> None:
        """Add new metabolites to a model.

        Will add a list of metabolites to the model object and add new
        constraints accordingly.

        The change is reverted upon exit when using the model as a context.

        Parameters
        ----------
        metabolite_list : list or Metabolite.
            A list of `cobra.core.Metabolite` objects. If it isn't an iterable
            container, the metabolite will be placed into a list.

        """
        if not hasattr(metabolite_list, "__iter__"):
            metabolite_list = [metabolite_list]
        if len(metabolite_list) == 0:
            return None

        # First check whether the metabolites exist in the model
        metabolite_list = [x for x in metabolite_list if x.id not in self.metabolites]

        bad_ids = [
            m for m in metabolite_list if not isinstance(m.id, str) or len(m.id) < 1
        ]
        if len(bad_ids) != 0:
            raise ValueError(f"invalid identifiers in {repr(bad_ids)}")

        for x in metabolite_list:
            x._model = self
        self.metabolites += metabolite_list

        # from cameo ...
        to_add = []
        for met in metabolite_list:
            if met.id not in self.constraints:
                constraint = self.problem.Constraint(Zero, name=met.id, lb=0, ub=0)
                to_add += [constraint]

        self.add_cons_vars(to_add)

        context = get_context(self)
        if context:
            context(partial(self.metabolites.__isub__, metabolite_list))
            for x in metabolite_list:
                # Do we care?
                context(partial(setattr, x, "_model", None))

    def remove_metabolites(
        self, metabolite_list: Union[List, Metabolite], destructive: bool = False
    ) -> None:
        """Remove a list of metabolites from the the object.

        The change is reverted upon exit when using the model as a context.

        Parameters
        ----------
        metabolite_list : list or Metaoblite
            A list of `cobra.core.Metabolite` objects. If it isn't an iterable
            container, the metabolite will be placed into a list.

        destructive : bool, optional
            If False then the metabolite is removed from all
            associated reactions.  If True then all associated
            reactions are removed from the Model (default False).
        """
        if not hasattr(metabolite_list, "__iter__"):
            metabolite_list = [metabolite_list]
        # Make sure metabolites exist in model
        metabolite_list = [x for x in metabolite_list if x.id in self.metabolites]
        for x in metabolite_list:
            x._model = None

            # remove reference to the metabolite in all groups
            associated_groups = self.get_associated_groups(x)
            for group in associated_groups:
                group.remove_members(x)

            if not destructive:
                for the_reaction in list(x._reaction):  # noqa W0212
                    the_coefficient = the_reaction._metabolites[x]  # noqa W0212
                    the_reaction.subtract_metabolites({x: the_coefficient})

            else:
                for x2 in list(x._reaction):  # noqa W0212
                    x2.remove_from_model()

        self.metabolites -= metabolite_list

        to_remove = [self.solver.constraints[m.id] for m in metabolite_list]
        self.remove_cons_vars(to_remove)

        context = get_context(self)
        if context:
            context(partial(self.metabolites.__iadd__, metabolite_list))
            for x in metabolite_list:
                context(partial(setattr, x, "_model", self))

    def add_boundary(
        self,
        metabolite: Metabolite,
        type: str = "exchange",
        reaction_id: Optional[str] = None,
        lb: Optional[float] = None,
        ub: Optional[float] = None,
        sbo_term: Optional[str] = None,
    ) -> Reaction:
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
        metabolites, i.e., a reversible reaction that adds or removes an
        intracellular metabolite.

        If you set the reaction `type` to something else, you must specify the
        desired identifier of the created reaction along with its upper and
        lower bound. The name will be given by the metabolite name and the
        given `type`.

        The change is reverted upon exit when using the model as a context.

        Parameters
        ----------
        metabolite : cobra.Metabolite
            Any given metabolite. The compartment is not checked but you are
            encouraged to stick to the definition of exchanges and sinks.
        type : {"exchange", "demand", "sink"}
            Using one of the pre-defined reaction types is easiest. If you
            want to create your own kind of boundary reaction choose
            any other string, e.g., 'my-boundary' (default "exchange").
        reaction_id : str, optional
            The ID of the resulting reaction. This takes precedence over the
            auto-generated identifiers but beware that it might make boundary
            reactions harder to identify afterwards when using `model.boundary`
            or specifically `model.exchanges` etc. (default None).
        lb : float, optional
            The lower bound of the resulting reaction (default None).
        ub : float, optional
            The upper bound of the resulting reaction (default None).
        sbo_term : str, optional
            A correct SBO term is set for the available types. If a custom
            type is chosen, a suitable SBO term should also be set (default None).

        Returns
        -------
        cobra.Reaction
            The created boundary reaction.

        Examples
        --------
        >>> from cobra.io load_model
        >>> model = load_model("textbook")
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
        ub = configuration.upper_bound if ub is None else ub
        lb = configuration.lower_bound if lb is None else lb
        types = {
            "exchange": ("EX", lb, ub, sbo_terms["exchange"]),
            "demand": ("DM", 0, ub, sbo_terms["demand"]),
            "sink": ("SK", lb, ub, sbo_terms["sink"]),
        }
        if type == "exchange":
            external = find_external_compartment(self)
            if metabolite.compartment != external:
                raise ValueError(
                    f"The metabolite is not an external metabolite (compartment is "
                    f"`{metabolite.compartment}` but should be `{external}`). "
                    f"Did you mean to add a demand or sink? If not, either change"
                    f" its compartment or rename the model compartments to fix this."
                )
        if type in types:
            prefix, lb, ub, default_term = types[type]
            if reaction_id is None:
                reaction_id = f"{prefix}_{metabolite.id}"
            if sbo_term is None:
                sbo_term = default_term
        if reaction_id is None:
            raise ValueError(
                "Custom types of boundary reactions require a custom "
                "identifier. Please set the `reaction_id`."
            )
        if reaction_id in self.reactions:
            raise ValueError(f"Boundary reaction '{reaction_id}' already exists.")
        name = f"{metabolite.name} {type}"
        rxn = Reaction(id=reaction_id, name=name, lower_bound=lb, upper_bound=ub)
        rxn.add_metabolites({metabolite: -1})
        if sbo_term:
            rxn.annotation["sbo"] = sbo_term
        self.add_reactions([rxn])
        return rxn

    def add_reactions(self, reaction_list: Iterable[Reaction]) -> None:
        """Add reactions to the model.

        Reactions with identifiers identical to a reaction already in the
        model are ignored.

        The change is reverted upon exit when using the model as a context.

        Parameters
        ----------
        reaction_list : list
            A list of `cobra.Reaction` objects
        """

        def existing_filter(rxn: Reaction) -> bool:
            """Check if the reaction does not exists in the model.

            Parameters
            ----------
            rxn: cobra.Reaction

            Returns
            -------
            bool
                False if reaction exists, True if it doesn't.
                If the reaction exists, will log a warning.
            """
            if rxn.id in self.reactions:
                logger.warning(f"Ignoring reaction '{rxn.id}' since it already exists.")
                return False
            return True

        # First check whether the reactions exist in the model.
        pruned = DictList(filter(existing_filter, reaction_list))

        context = get_context(self)

        # Add reactions. Also take care of genes and metabolites in the loop.
        for reaction in pruned:
            reaction._model = self
            if context:
                context(partial(setattr, reaction, "_model", None))
            # Build a `list()` because the dict will be modified in the loop.
            for metabolite in list(reaction.metabolites):
                # TODO: Maybe this can happen with
                #  Reaction.add_metabolites(combine=False)
                # TODO: Should we add a copy of the metabolite instead?
                if metabolite not in self.metabolites:
                    self.add_metabolites(metabolite)
                # A copy of the metabolite exists in the model, the reaction
                # needs to point to the metabolite in the model.
                else:
                    # FIXME: Modifying 'private' attributes is horrible.
                    stoichiometry = reaction._metabolites.pop(metabolite)
                    model_metabolite = self.metabolites.get_by_id(metabolite.id)
                    reaction._metabolites[model_metabolite] = stoichiometry
                    model_metabolite._reaction.add(reaction)
                    if context:
                        context(partial(model_metabolite._reaction.remove, reaction))
            reaction.update_genes_from_gpr()

        self.reactions += pruned

        if context:
            context(partial(self.reactions.__isub__, pruned))

        # from cameo ...
        self._populate_solver(pruned)

    def remove_reactions(
        self,
        reactions: Union[str, Reaction, List[Union[str, Reaction]]],
        remove_orphans: bool = False,
    ) -> None:
        """Remove reactions from the model.

        The change is reverted upon exit when using the model as a context.

        Parameters
        ----------
        reactions : list or reaction or str
            A list with reactions (`cobra.Reaction`), or their id's, to remove.
            Reaction will be placed in a list. Str will be placed in a list and used to
            find the reaction in the model.
        remove_orphans : bool, optional
            Remove orphaned genes and metabolites from the model as
            well (default False).
        """
        if isinstance(reactions, str) or hasattr(reactions, "id"):
            warn("need to pass in a list")
            reactions = [reactions]

        context = get_context(self)

        for reaction in reactions:

            # Make sure the reaction is in the model
            try:
                reaction = self.reactions[self.reactions.index(reaction)]
            except ValueError:
                warn(f"{reaction} not in {self}")

            else:
                forward = reaction.forward_variable
                reverse = reaction.reverse_variable

                if context:

                    obj_coef = reaction.objective_coefficient

                    if obj_coef != 0:
                        context(
                            partial(
                                self.solver.objective.set_linear_coefficients,
                                {forward: obj_coef, reverse: -obj_coef},
                            )
                        )

                    context(partial(self._populate_solver, [reaction]))
                    context(partial(setattr, reaction, "_model", self))
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

    def add_groups(self, group_list: Union[str, Group, List[Group]]) -> None:
        """Add groups to the model.

        Groups with identifiers identical to a group already in the model are
        ignored.

        If any group contains members that are not in the model, these members
        are added to the model as well. Only metabolites, reactions, and genes
        can have groups.

        Parameters
        ----------
        group_list : list or str or Group
            A list of `cobra.Group` objects to add to the model. Can also be a single
            group or a string representing group id. If the input is not a list, a
            warning is raised.
        """

        def existing_filter(new_group: Group) -> bool:
            """Check if the group does not exist.

            Parameters
            ----------
            new_group: cobra.Group
                Group to check.

            Returns
            -------
            bool
                False if the group already exists, True if it doesn't.
            """
            if new_group.id in self.groups:
                logger.warning(
                    f"Ignoring group '{new_group.id}'" f" since it already exists."
                )
                return False
            return True

        if isinstance(group_list, str) or hasattr(group_list, "id"):
            warn("need to pass in a list")
            group_list = [group_list]

        pruned = DictList(filter(existing_filter, group_list))

        for group in pruned:
            group._model = self
            for member in group.members:
                # If the member is not associated with the model, add it
                if isinstance(member, Metabolite) and member not in self.metabolites:
                    self.add_metabolites([member])
                if isinstance(member, Reaction) and member not in self.reactions:
                    self.add_reactions([member])
                # TODO(midnighter): `add_genes` method does not exist.
                # if isinstance(member, Gene):
                #     if member not in self.genes:
                #         self.add_genes([member])

            self.groups += [group]

    def remove_groups(self, group_list: Union[str, Group, List[Group]]) -> None:
        """Remove groups from the model.

        Members of each group are not removed
        from the model (i.e. metabolites, reactions, and genes in the group
        stay in the model after any groups containing them are removed).

        Parameters
        ----------
        group_list : list or str or Group
            A list of `cobra.Group` objects to remove from the model. Can also be a
            single group or a string representing group id. If the input is not a list,
             a warning is raised.
        """
        if isinstance(group_list, str) or hasattr(group_list, "id"):
            warn("need to pass in a list")
            group_list = [group_list]

        for group in group_list:
            # make sure the group is in the model
            if group.id not in self.groups:
                logger.warning(f"{group!r} not in {self!r}. Ignored.")
            else:
                self.groups.remove(group)
                group._model = None

    def get_associated_groups(
        self, element: Union[Reaction, Gene, Metabolite]
    ) -> List[Group]:
        """Get list of groups for element.

        Returns a list of groups that an element (reaction, metabolite, gene)
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

    def add_cons_vars(
        self, what: Union[List["CONS_VARS"], Tuple["CONS_VARS"]], **kwargs
    ) -> None:
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

    def remove_cons_vars(
        self, what: Union[List["CONS_VARS"], Tuple["CONS_VARS"]]
    ) -> None:
        """Remove variables and constraints from problem.

        Remove variables and constraints from the model's mathematical
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
    def problem(self) -> "optlang.interface":
        """Get the interface to the model's underlying mathematical problem.

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
    def variables(self) -> "Container":
        """Get the mathematical variables in the cobra model.

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
    def constraints(self) -> "Container":
        """Get the constraints in the cobra model.

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
    def boundary(self) -> List[Reaction]:
        """Boundary reactions in the model.

        Reactions that either have no substrate or product.

        Returns
        -------
        list
            A list of reactions that either have no substrate or product and
            only one metabolite overall.
        """
        return [rxn for rxn in self.reactions if rxn.boundary]

    @property
    def exchanges(self) -> List[Reaction]:
        """Exchange reactions in model.

        Reactions that exchange mass with the exterior. Uses annotations
        and heuristics to exclude non-exchanges such as sink reactions.

        Returns
        -------
        list
            A list of reactions that satisfy the conditions for exchange reactions.

        See Also
        --------
        cobra.medium.find_boundary_types
        """
        return find_boundary_types(self, "exchange", None)

    @property
    def demands(self) -> List[Reaction]:
        """Demand reactions in model.

        Irreversible reactions that accumulate or consume a metabolite in
        the inside of the model.

        Returns
        -------
        list
            A list of reactions that are demand reactions (reactions that
            accumulate/consume a metabolite irreversibly).

        See Also
        --------
        cobra.medium.find_boundary_types
        """
        return find_boundary_types(self, "demand", None)

    @property
    def sinks(self) -> List[Reaction]:
        """Sink reactions in model.

        Reversible reactions that accumulate or consume a metabolite in
        the inside of the model.

        Returns
        -------
        list
            A list of reactions that are demand reactions (reactions that
            accumulate/consume a metabolite reversibly).

        See Also
        --------
        cobra.medium.find_boundary_types
        """
        return find_boundary_types(self, "sink", None)

    def _populate_solver(
        self,
        reaction_list: List[Reaction],
        metabolite_list: Optional[List[Metabolite]] = None,
    ) -> None:
        """Populate attached solver with constraints and variables.

        Populate attached solver with constraints and variables that
        model the provided reactions.

        Parameters
        ----------
        reaction_list: list
            A list of cobra.Reaction to add to the solver. This list will be
            constrained.
        metabolite_list: list, optional
            A list of cobra.Metabolite  to add to the solver. This list will be
            constrained (default None).
        """
        constraint_terms = AutoVivification()
        to_add = []
        if metabolite_list is not None:
            for met in metabolite_list:
                to_add += [self.problem.Constraint(Zero, name=met.id, lb=0, ub=0)]
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
            for metabolite, coeff in reaction.metabolites.items():
                if metabolite.id in self.constraints:
                    constraint = self.constraints[metabolite.id]
                else:
                    constraint = self.problem.Constraint(
                        Zero, name=metabolite.id, lb=0, ub=0
                    )
                    self.add_cons_vars(constraint, sloppy=True)
                constraint_terms[constraint][forward_variable] = coeff
                constraint_terms[constraint][reverse_variable] = -coeff

        self.solver.update()
        for reaction in reaction_list:
            reaction = self.reactions.get_by_id(reaction.id)
            reaction.update_variable_bounds()
        for constraint, terms in constraint_terms.items():
            constraint.set_linear_coefficients(terms)

    def slim_optimize(
        self, error_value: Optional[float] = float("nan"), message: Optional[str] = None
    ) -> float:
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
           optimization fails (default float("nan")).
        message : str, optional
           Error message to use if the model optimization did not succeed (default
           None).

        Returns
        -------
        float
            The objective value. Returns the error value if optimization failed and
            error_value was not None.

        Raises
        ------
        OptimizationError
            If error_value was set as None and the optimization fails.
        """
        self.solver.optimize()
        if self.solver.status == optlang.interface.OPTIMAL:
            return self.solver.objective.value
        elif error_value is not None:
            return error_value
        else:
            assert_optimal(self, message)

    def optimize(
        self, objective_sense: Optional[str] = None, raise_error: bool = False
    ) -> "Solution":
        """
        Optimize the model using flux balance analysis.

        Parameters
        ----------
        objective_sense : {None, 'maximize' 'minimize'}, optional
            Whether fluxes should be maximized or minimized. In case of None,
            the previous direction is used (default None).
        raise_error : bool
            If true, raise an OptimizationError if solver status is not
             optimal (default False).

        Returns
        -------
        Solution

        Notes
        -----
        Only the most commonly used parameters are presented here.  Additional
        parameters for cobra.solvers may be available and specified with the
        appropriate keyword argument.

        """
        original_direction = self.objective.direction
        self.objective.direction = {"maximize": "max", "minimize": "min"}.get(
            objective_sense, original_direction
        )
        self.slim_optimize()
        solution = get_solution(self, raise_error=raise_error)
        self.objective.direction = original_direction
        return solution

    def repair(
        self, rebuild_index: bool = True, rebuild_relationships: bool = True
    ) -> None:
        """Update all indexes and pointers in a model.

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
                rxn.update_genes_from_gpr()
                for met in rxn._metabolites:
                    met._reaction.add(rxn)

        # point _model to self
        for dict_list in (self.reactions, self.genes, self.metabolites, self.groups):
            for entity in dict_list:
                entity._model = self

    @property
    def objective(self) -> Union[optlang.Objective]:
        """Get the solver objective.

        With optlang, the objective is not limited to a simple linear summation of
        individual reaction fluxes, making the return value ambiguous.

        Henceforth, use `cobra.util.solver.linear_reaction_coefficients` to
        get a dictionary of reactions with their linear coefficients (empty
        if there are none).

        """
        return self.solver.objective

    @objective.setter
    def objective(
        self,
        value: Union[
            Dict[Reaction, "Container"],
            str,
            int,
            Reaction,
            "optlang.interface.Objective",
            Basic,
        ],
    ) -> None:
        """Set the solver objective.

        Parameters
        ----------
        value: dict or str or int or Reaction or  optlang.interface.Container, Reaction
        or Basic
            The set value can be dictionary (reactions as keys, linear coefficients as
            values), string (reaction identifier), int (reaction index), Reaction or
            problem.Objective or sympy expression directly interpreted as objectives.

        When using in a context, this attribute can be set temporarily.
        """
        if isinstance(value, Basic):
            value = self.problem.Objective(value, sloppy=False)
        if not isinstance(value, (dict, optlang.interface.Objective)):
            try:
                reactions = self.reactions.get_by_any(value)
            except KeyError:
                raise ValueError("invalid objective")
            value = {rxn: 1 for rxn in reactions}
            # TODO - check that it is reset with context.
        set_objective(self, value, additive=False)

    @property
    def objective_direction(self) -> str:
        """Get the objective direction.

        Returns
        -------
        str
            Objective direction as string. Should be "max" or "min".
        """
        return self.solver.objective.direction

    @objective_direction.setter
    @resettable
    def objective_direction(self, value: str) -> None:
        """Set the objective direction.

        When used in a context, this attribute is set temporarily.

        Parameters
        ----------
        value: {"max", "min"}
            String of "max" or "min" for direction.

        Raises
        ------
        ValueError
            If given direction isn't max or min.
        """
        value = value.lower()
        if value.startswith("max"):
            self.solver.objective.direction = "max"
        elif value.startswith("min"):
            self.solver.objective.direction = "min"
        else:
            raise ValueError(f"Unknown objective direction '{value}'.")

    def summary(
        self,
        solution: Optional["Solution"] = None,
        fva: Union["pd.DataFrame", float, None] = None,
    ) -> "ModelSummary":
        """
        Create a summary of the exchange fluxes of the model.

        Parameters
        ----------
        solution : cobra.Solution, optional
            A previous model solution to use for generating the summary. If
            ``None``, the summary method will generate a parsimonious flux
            distribution (default None).
        fva : pd.DataFrame or float, optional
            Whether or not to include flux variability analysis in the output.
            If given, `fva` should either be a previous FVA solution matching the
            model or a float between 0 and 1 representing the fraction of the
            optimum objective to be searched (default None).

        Returns
        -------
        cobra.ModelSummary

        See Also
        --------
        Reaction.summary
        Metabolite.summary

        """
        from cobra.summary import ModelSummary

        return ModelSummary(model=self, solution=solution, fva=fva)

    def __enter__(self) -> "Model":
        """Record future changes to the model.

        Record all future changes to the model, undoing them when a call to
        __exit__ is received. Creates a new context and adds it to the stack.

        Returns
        -------
        cobra.Model
            Returns the model with context added.
        """
        try:
            self._contexts.append(HistoryManager())
        except AttributeError:
            self._contexts = [HistoryManager()]

        return self

    def __exit__(self, type, value, traceback) -> None:
        """Pop the top context manager and trigger the undo functions."""
        context = self._contexts.pop()
        context.reset()

    def merge(
        self,
        right: "Model",
        prefix_existing: Optional[str] = None,
        inplace: bool = True,
        objective: str = "left",
    ) -> "Model":
        """Merge two models to create a model with the reactions from both models.

        Custom constraints and variables from right models are also copied
        to left model, however note that, constraints and variables are
        assumed to be the same if they have the same name.

        Parameters
        ----------
        right : cobra.Model
            The model to add reactions from
        prefix_existing : string or optional
            Prefix the reaction identifier in the right that already exist
            in the left model with this string (default None).
        inplace : bool
            Add reactions from right directly to left model object.
            Otherwise, create a new model leaving the left model untouched.
            When done within the model as context, changes to the models are
            reverted upon exit (default True).
        objective : {"left", "right", "sum"}
            One of "left", "right" or "sum" for setting the objective of the
            resulting model to that of the corresponding model or the sum of
            both (default "left").

        Returns
        -------
        cobra.Model
            The merged model.
        """
        if inplace:
            new_model = self
        else:
            new_model = self.copy()
            new_model.id = f"{self.id}_{right.id}"
        new_reactions = deepcopy(right.reactions)
        if prefix_existing is not None:
            existing = new_reactions.query(lambda rxn: rxn.id in self.reactions)
            for reaction in existing:
                reaction.id = f"{prefix_existing}{reaction.id}"
        new_model.add_reactions(new_reactions)
        interface = new_model.problem
        new_vars = [
            interface.Variable.clone(v)
            for v in right.variables
            if v.name not in new_model.variables
        ]
        new_model.add_cons_vars(new_vars)
        new_cons = [
            interface.Constraint.clone(c, model=new_model.solver)
            for c in right.constraints
            if c.name not in new_model.constraints
        ]
        new_model.add_cons_vars(new_cons, sloppy=True)
        new_model.objective = dict(
            left=self.objective,
            right=right.objective,
            sum=self.objective.expression + right.objective.expression,
        )[objective]
        return new_model

    def _repr_html_(self) -> str:
        """Get HTML represenation of the model.

        Returns
        -------
        str
            Model representation as HTML string.
        """
        return f"""
        <table>
            <tr>
                <td><strong>Name</strong></td>
                <td>{self.id}</td>
            </tr><tr>
                <td><strong>Memory address</strong></td>
                <td>{f"{id(self):x}"}</td>
            </tr><tr>
                <td><strong>Number of metabolites</strong></td>
                <td>{len(self.metabolites)}</td>
            </tr><tr>
                <td><strong>Number of reactions</strong></td>
                <td>{len(self.reactions)}</td>
            </tr><tr>
                <td><strong>Number of genes</strong></td>
                <td>{len(self.genes)}</td>
            </tr><tr>
                <td><strong>Number of groups</strong></td>
                <td>{len(self.groups)}</td>
            </tr><tr>
                <td><strong>Objective expression</strong></td>
                <td>{format_long_string(str(self.objective.expression), 100)}</td>
            </tr><tr>
                <td><strong>Compartments</strong></td>
                <td>{", ".join(v if v else k for k, v in
                               self.compartments.items())}</td>
            </tr>
          </table>"""
