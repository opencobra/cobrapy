# -*- coding: utf-8 -*-

from __future__ import absolute_import

from six import string_types


class Object(object):
    """Defines common behavior of object in cobra.core"""

    def __init__(self, id=None, name=""):
        """A simple object with an identifier

        Parameters
        ----------
        id: None or a string
            the identifier to associate with the object
        """
        self._id = id
        self.name = name

        self.notes = {}
        self.annotation = {}

    @property
    def id(self):
        return getattr(self, "_id", None)

    @id.setter
    def id(self, value):
        if value == self.id:
            pass
        elif not isinstance(value, string_types):
            raise TypeError("ID must be a string")
        elif getattr(self, "_model", None) is not None:
            self._set_id_with_model(value)
        else:
            self._id = value

    def _set_id_with_model(self, value):
        self._id = value

    def __getstate__(self):
        """To prevent excessive replication during deepcopy."""
        state = self.__dict__.copy()
        if '_model' in state:
            state['_model'] = None
        return state

    def __repr__(self):
        return "<%s %s at 0x%x>" % (self.__class__.__name__, self.id, id(self))

    def __str__(self):
        return str(self.id)
