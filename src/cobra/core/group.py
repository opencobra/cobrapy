"""Define the group class."""

from typing import Iterable, Optional, Set
from warnings import warn

from .object import Object


class Group(Object):
    """
    Manage groups via this implementation of the SBML group specification.

    `Group` is a class for holding information regarding a pathways,
    subsystems, or other custom groupings of objects within a cobra.Model
    object.

    Parameters
    ----------
    id : str
        The identifier to associate with this group
    name : str, optional
        A human readable name for the group
    members : iterable, optional
        A DictList containing references to cobra.Model-associated objects
        that belong to the group.
    kind : {"collection", "classification", "partonomy"}, optional
        The kind of group, as specified for the Groups feature in the SBML
        level 3 package specification. Can be any of "classification",
        "partonomy", or "collection". The default is "collection".
        Please consult the SBML level 3 package specification to ensure you
        are using the proper value for kind. In short, members of a
        "classification" group should have an "is-a" relationship to the group
        (e.g. member is-a polar compound, or member is-a transporter).
        Members of a "partonomy" group should have a "part-of" relationship
        (e.g. member is part-of glycolysis). Members of a "collection" group
        do not have an implied relationship between the members, so use this
        value for kind when in doubt (e.g. member is a gap-filled reaction,
        or member is involved in a disease phenotype).
    """

    KIND_TYPES = ("collection", "classification", "partonomy")

    def __init__(
        self,
        id: str,
        name: str = "",
        members: Optional[Iterable] = None,
        kind: Optional[str] = None,
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

        self._members = set() if members is None else set(members)
        self._kind = None
        self.kind = "collection" if kind is None else kind
        # self.model is None or refers to the cobra.Model that
        # contains self
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
    def members(self) -> Set:
        """Get members of the group.

        Returns
        -------
        Set
            A Set containing the members of the group.
        """
        return self._members

    @property
    def kind(self) -> str:
        """Return the group kind.

        Returns
        -------
        str
            The group kind. Should be one of the three types allowed in SBML.

        """
        return self._kind

    @kind.setter
    def kind(self, kind: str) -> None:
        """Set the group kind.

        Parameters
        ----------
        kind: str
            Must be one of the allowed kind types "collection", "classification",
            "partonomy".
            If kind is not one of these types, a ValueError Exception is raised.
        """
        kind = kind.lower()
        if kind in self.KIND_TYPES:
            self._kind = kind
        else:
            raise ValueError(f"Kind can only by one of: {', '.join(self.KIND_TYPES)}.")

    def add_members(self, new_members: list) -> None:
        """Add objects to the group.

        Parameters
        ----------
        new_members : list
            A list of cobrapy objects to add to the group.
        """
        if isinstance(new_members, str) or hasattr(new_members, "id"):
            warn("need to pass in a list")
            new_members = [new_members]

        self._members.update(new_members)

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

        self._members.difference_update(to_remove)
