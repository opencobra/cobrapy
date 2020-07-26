# -*- coding: utf-8 -*-

from __future__ import absolute_import

from collections import defaultdict

from six import string_types

from cobra.core.metadata import CVList, MetaData, Notes


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

        self._notes = Notes()
        self._annotation = MetaData()

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

    @property
    def notes(self):
        return getattr(self, "_notes", None)

    @notes.setter
    def notes(self, data):
        if isinstance(data, Notes):
            self._notes = data
            return
        self._notes.notes_xhtml = data

    @property
    def annotation(self):
        return getattr(self, "_annotation", None)

    @annotation.setter
    def annotation(self, value):
        if not (isinstance(value, dict) or isinstance(value, MetaData)):
            raise TypeError("The data passed for annotation must be inside "
                            "a dictionary: {}".format(value))
        else:
            if isinstance(value, MetaData):
                self._annotation = value
            else:
                self._annotation.cvterms._annotations = defaultdict(list)
                self._annotation.cvterms._cvterms = defaultdict(CVList)
                self._annotation.cvterms.add_simple_annotations(value)

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
