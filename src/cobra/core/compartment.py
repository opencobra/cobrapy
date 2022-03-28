"""Define the group class."""

from typing import Iterable, Optional, Union, List, FrozenSet
from warnings import warn

from .dictlist import DictList
from .group import Group
from .. import Metabolite, Reaction, Model


class Compartment(Group):
    """
    Manage groups via this implementation of the SBML group specification.

    `Compartment` is a class for holding information regarding a bounded space in
    which species are located. You can think of it as a location in a model,
    usually representing organelles like the Endoplasmic Reticulum (ER).

    Parameters
    ----------
    id : str
        The identifier to associate with this group
    name : str, optional
        A human readable name for the group
    members : iterable, optional
        A DictList containing references to cobra.Model-associated objects
        that belong to the group. Members should be metabolites or genes, where
        reactions are inferred from metabolites.
    dimensions: float, optional
        Compartments can have dimensions defined, if they are volume (3 dimensions) or
        2 (a two-dimensional compartment, a surface, like a membrane). In theory, this
        can be 1 or 0 dimensions, and even partial dimensions, but that will needlessly
        complicate the math. The number of dimensions influences the size and units
        used.
    """
    def __init__(self, id: str, name: str = "", members: Optional[Iterable] = None,
                 dimensions: Optional[float] = None):
        """Initialize the group object.

         id : str
            The identifier to associate with this group
        name : str, optional
            A human readable name for the group
        members : iterable, optional
            A DictList containing references to cobra.Model-associated objects
            that belong to the group.
        dimensions: float, optional
            Compartments can have dimensions defined, if they are volume (3 dimensions) or
            2 (a two-dimensional compartment, a surface, like a membrane). In theory, this
            can be 1 or 0 dimensions, and even partial dimensions. The number of
            dimensions influences the size and units used.

        Raises
        ------
        TypeError if given anything other than metaoblites or reactions.
        """
        super().__init__(id, name, members)
        for x in members:
            if not isinstance(x, (Metabolite, Reaction)):
                raise(TypeError, f"Compartments should have only "
                                 f"reactions or metabolites. {x} is a {type(x)}.")
        self._members = DictList() if members is None else DictList(members)
        self._compartment_type = None
        self.__delattr__("kind") # Compartments don't have kind
        # self.model is None or refers to the cobra.Model that
        # contains self
        self._dimensions = dimensions
        self._model = None

    def add_members(self, new_members: Union[str, Metabolite, Reaction,
                                             List[Union[Reaction, Metabolite]]]) -> None:
        """Add objects to the compartment.

        Parameters
        ----------
        new_members : list or str or Metabolite or Reaction
            A list of cobrapy Metabolites or Genes to add to the group. If it isn't a
            list a warning will be raised.
            If it isn't a metaoblite or a gene, an error will be raised.

        Raises
        ------
        TypeError - if given any object other than Metaoblite or Reaction
        """
        if isinstance(new_members, str) or hasattr(new_members, "id"):
            warn("need to pass in a list")
            new_members = [new_members]

        #TODO - add some filtering on type. What happens if given a string? Check
        # DictList and groups.

        self._members.union(new_members)
        for _member in new_members:
            _member._compartment = self

    def remove_members(self, to_remove: list) -> None:
        """Remove objects from the group.

        Parameters
        ----------
        to_remove : list
            A list of cobra objects to remove from the group
        """
        if isinstance(to_remove, str) or hasattr(to_remove, "id"):
            warn("need to pass in a list")
            to_remove = [to_remove]

        for member_to_remove in to_remove:
            self._members.remove(member_to_remove)
            member_to_remove._compartment = None

    @property
    def metabolites(self) -> DictList[Metabolite]:
        """Return the metaoblites in the compartment, if any.

        Returns
        -------
        DictList:
            DictList of metaoblites if any are present in the compartment. If there are
            no metaoblties, will return an empty DictList.

        """
        return self._members.query(lambda x: isinstance(x, Metabolite))

    @property
    def inferred_reactions(self) -> FrozenSet[Reaction]:
        """Return the reactions whose metabolites are in the compartment.

        This is returned as a FrozenSet of reactions for each metabolite in the
        compartment, if any.

        Returns
        -------
        FrozenSet of cobra.Reactions
            Reactions that have metabolites that belong to this compartment.

        Returns
        -------
        FrozenSet of cobra.Reactions
            Reactions that have metabolites that belong to this compartment.
        """
        rxn_set = set()
        for met in self.metabolites:
            rxn_set.update(met._reactions)
        return frozenset(rxn_set)

    @property
    def assigned_reactions(self) -> FrozenSet[Reaction]:
        """Return the reactions who were assigned to this compartment.

        This is returned as a FrozenSet of reactions for each metabolite in the
        compartment, if any.

        Returns
        -------
        FrozenSet of cobra.Reactions
            Reactions that have metabolites that belong to this compartment.

        Returns
        -------
        FrozenSet of cobra.Reactions
            Reactions that were assigned to this compartment, if any.
        """
        return frozenset(self._members.query(lambda x: isinstance(x, Reaction)))

    @property
    def reactions(self) -> FrozenSet[Reaction]:
        """Return the reactions who belong to this compartment.

        This is returned as a FrozenSet of reactions for each metabolite in the
        compartment, if any, and the reactions that were assigned to this compartment
        directly.

        Returns
        -------
        FrozenSet of cobra.Reactions
            Reactions that belong to this compartment, both assigned and inferred.
        """
        direct_set = set(self.assigned_reactions)
        return frozenset(direct_set.union(self.inferred_reactions))

    def __contains__(self, member: Union[Metabolite, Reaction]):
        return self.members.__contains__(member)

    def merge(self, other):
        warn("Not implemented yet")
        return

    def copy(self, new_id: str, new_model: Model, new_name: Optional[str] = None):
        warn("Not implemented yet")
        return
