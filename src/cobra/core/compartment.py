"""Define the group class."""

from typing import Iterable, Optional, Union, List, FrozenSet
from warnings import warn

from .dictlist import DictList
from .object import Object
from .. import Metabolite, Gene, Reaction, Model


class Compartment(Object):
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
    compartment_type : str, optional
        SBML Level 2 Versions 2–4 provide the compartment type as a grouping construct
        that can be used to establish a relationship between multiple Compartment
        objects. A CompartmentType object only has an identity, and this identity can
        only be used to indicate that particular Compartment objects in the model
        belong to this type. This may be useful for conveying a modeling intention,
        such as when a model contains many similar compartments, either by their
        biological function or the reactions they carry. Without a compartment type
        construct, it would be impossible within SBML itself to indicate that all of
        the compartments share an underlying conceptual relationship because each
        SBML compartment must be given a unique and separate identity. A
        CompartmentType has no mathematical meaning in SBML—it has no effect on
        a model's mathematical interpretation. Simulators and other numerical analysis
        software may ignore CompartmentType definitions and references in a model.
    dimensions: float, optional
        Compartments can have dimensions defined, if they are volume (3 dimensions) or
        2 (a two-dimensional compartment, a surface, like a membrane). In theory, this
        can be 1 or 0 dimensions, and even partial dimensions, but that will needlessly
        complicate the math. The number of dimensions influences the size and units
        used.
    """
    def __init__(
        self,
        id: str,
        name: str = "",
        members: Optional[Iterable] = None,
        compartment_type: Optional[str] = None,
        dimensions: Optional[float] = None,
    ):
        """Initialize the group object.

         id : str
            The identifier to associate with this group
        name : str, optional
            A human readable name for the group
        members : iterable, optional
            A DictList containing references to cobra.Model-associated objects
            that belong to the group.
        kind : {"collection", "classification", "partonomy"}, optional
            The kind of group, as specified for the Groups feature in the SBML
            level 3 package specification.
        """
        Object.__init__(self, id, name)

        self._members = DictList() if members is None else DictList(members)
        self._compartment_type = None
        self.compartment_type = "" if compartment_type is None else compartment_type
        # self.model is None or refers to the cobra.Model that
        # contains self
        self._dimensions = dimensions
        self._model = None

    def __len__(self) -> int:
        """Get length of group.

        Returns
        -------
        int
            An int with the length of the group.

        """
        return len(self._members)

    # read-only
    @property
    def members(self) -> DictList:
        """Get members of the group.

        Returns
        -------
        DictList
            A dictlist containing the members of the group.
        """
        return self._members

    @property
    def compartment_type(self) -> str:
        """Return the compartment type.

        Returns
        -------
        str
            The compartment type.

        """
        return self._compartment_type

    @compartment_type.setter
    def compartment_type(self, compartment_type: str) -> None:
        """Set the compartment type.

        Parameters
        ----------
        compartment_type: str

        """
        self._compartment_type = compartment_type

    def add_members(self, new_members: Union[str, Metabolite, Reaction,
                                             List[Union[Reaction, Metabolite]]]) -> None:
        """Add objects to the group.

        Parameters
        ----------
        new_members : list or str or Metabolite or Reaction
            A list of cobrapy Metabolites or Genes to add to the group. If it isn't a
            list a warning will be raised.
            If it isn't a metaoblite or a gene, an error will be raised.

        Raises
        ------
        TypeError - if given any object other than Metaoblite or Gene
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
    def reactions(self) -> Optional[FrozenSet[Reaction]]:
        """Return the reactions whose metabolites are in the compartment.

        This is returned as a FrozenSet of reactions for each metabolite in the
        compartment, if any.

        Returns
        -------
        FrozenSet of cobra.Reactions
            Reactions that have metabolites that belong to this compartment.
        """
        direct_set = set(self._members.query(lambda x: isinstance(x, Reaction)))
        rxn_set = set()
        for met in self.metabolites:
            rxn_set.update(met._reactions)
        return frozenset(rxn_set.union(direct_set))

    def __contains__(self, member: Union[Metabolite, Gene]):
        return member.compartment is self

    def merge(self, other):
        warn("Not implemented yet")
        return

    def copy(self, new_id: str, new_model: Model, new_name: Optional[str] = None):
        warn("Not implemented yet")
        return
