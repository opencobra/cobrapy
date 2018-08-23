# -*- coding: utf-8 -*-

from __future__ import absolute_import

from cobra.core.object import Object
from six import string_types
from warnings import warn

kind_types = ["collection", "classification", "partonomy"]


class Group(Object):
    """Group is a class for holding information regarding
    a pathways, subsystems, or other custom groupings of objects
    within a cobra.Model object.

    Parameters
    ----------
    id : string
        The identifier to associate with this group
    name : string
        A human readable name for the group
    members : list
        A list object containing references to cobra.Model-associated objects
        that belong to the group.
    kind : string
        The kind of group, as specified for the Groups feature in the SBML
        level 3 package specification. Can be any of "classification",
        "partonomy", or "collection". Please consult the SBML level 3 package
        specification to ensure you are using the proper value for kind. In
        short, members of a "classification" group should have an "is-a"
        relationship to the group (e.g. member is-a polar compound, or
        member is-a transporter). Members of a "partonomy" group should have a
        "part-of" relationship (e.g. member is part-of glycolysis). Members of
        a "collection" group do not have an implied relationship between the
        members, so use this value for kind when in doubt (e.g. member is a
        gap-filled reaction, or member is involved in a disease phenotype).
    """

    def __init__(self, id=None, name='', members=[], kind=''):
        Object.__init__(self, id, name)

        self._members = set(members)
        self._kind = kind

        # self.model is None or refers to the cobra.Model that
        # contains self
        self._model = None

    # read-only
    @property
    def members(self):
        return getattr(self, "_members", None)

    @members.setter
    def members(self, members):
        self._members = set(members)

    @property
    def kind(self):
        return getattr(self, "_kind", '')

    @kind.setter
    def kind(self, kind):
        if kind in kind_types:
            self._kind = kind
        else:
            raise ValueError("kind can only by one of: " + str(kind_types))

    @property
    def model(self):
        """returns the model the group is a part of"""
        return self._model

    def add_members(self, members_list):
        """Add objects to the group.

        Parameters
        ----------
        members_to_add : list
            list of cobrapy objects to add to the group.
        """

        if isinstance(members_list, string_types) or \
                hasattr(members_list, "id"):
            warn("need to pass in a list")
            members_list = [members_list]

        new_members = []
        _id_to_members = dict([(x.id, x) for x in self._members])

        # Check for duplicate members in the group
        for member in members_list:
            # we only need to add the member if it ins't already in the group
            if member.id not in _id_to_members:
                new_members.append(member)

        self._members = self._members.union(set(new_members))

    def remove(self, members_list):
        """Remove objects from the group.

        Parameters
        ----------
        members_to_remove : list
            list of cobrapy objects to remove from the group
        """

        if isinstance(members_list, string_types) or \
                hasattr(members_list, "id"):
            warn("need to pass in a list")
            members_list = [members_list]

        for member in members_list:
            self._members.discard(member)
