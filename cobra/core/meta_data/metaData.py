# -*- coding: utf-8 -*-

from __future__ import absolute_import

from cobra.core.meta_data.cvTerm import CVTerm
from cobra.core.meta_data.history import History
from cobra.core.meta_data.keyValuePair import ListOfKeyValue


class MetaData:
    """Class representation of the meta-data of the component.
    It is a combination of three classes i.e CVTerm, History
    and KeyValuePair class.

    Parameters
    ----------
    cvterm : dict, CVTerm
        The cvterm holds data for external resources
    history : dict, History
        The history is holding the data about the creator,
        created and modified dates.
    listofKeyValue : list
        Some key-value pairs which are not suitable to be
        represented anywhere else in the model.

    """

    def __init__(self, cvterm=None, history=None, listofKeyValue=None):
        # setting the cvterm
        if cvterm is None:
            self._cvTerms = CVTerm()
        elif isinstance(cvterm, CVTerm):
            self._cvTerms = cvterm
        elif isinstance(cvterm, dict):
            self._cvTerms = CVTerm(cvterm)
        else:
            raise TypeError("Invalid format passed for cvterm")
        # setting the history of the component
        if history is None:
            self._history = History()
        elif isinstance(history, History):
            self._history = history
        elif isinstance(history, dict):
            if "creator" not in history:
                history["creator"] = []
            if "created" not in history:
                history["created"] = None
            if "modified" not in history:
                history["modified"] = []
            self._history = History(history["creator"],
                                    history["created"], history["modified"])
        else:
            raise TypeError("Invalid format passed for history")
        # setting the list of key-value pair
        if listofKeyValue is not None:
            if isinstance(listofKeyValue, ListOfKeyValue):
                self._listofKeyValue = listofKeyValue
            elif isinstance(listofKeyValue, list):
                self._listofKeyValue = ListOfKeyValue(listofKeyValue)
            else:
                raise TypeError("Key value pairs must be passed in a list")
        else:
            self._listofKeyValue = ListOfKeyValue()

    @property
    def cvTerms(self):
        return getattr(self, "_cvTerms", None)

    @cvTerms.setter
    def cvTerms(self, value):
        if value == self._cvTerms:
            pass
        elif isinstance(value, CVTerm):
            self._cvTerms = value
        elif isinstance(value, dict):
            self._cvTerms = CVTerm(cvterm)
        else:
            raise TypeError("This passed format for cvterm is invalid")

    @property
    def history(self):
        return getattr(self, "_history", None)

    @history.setter
    def history(self, value):
        if value == self._history:
            pass
        elif isinstance(value, History):
            self._history = value
        elif isinstance(value, dict):
            if "creator" not in history:
                history["creator"] = []
            if "created" not in history:
                history["created"] = None
            if "modified" not in history:
                history["modified"] = []
            self._history = History(value["creator"],
                                    value["created"], value["modified"])
        else:
            raise TypeError("This passed format for history is invalid")

    @property
    def listofKeyValue(self):
        return getattr(self, "_listofKeyValue", [])

    @listofKeyValue.setter
    def listofKeyValue(self, value):
        if value == self._listofKeyValue:
            pass
        elif isinstance(value, ListOfKeyValue):
            self._listofKeyValue = value
        elif isinstance(value, list):
            self._listofKeyValue = ListOfKeyValue(value)
        else:
            raise TypeError("This passed format for listofKeyValue is "
                            "invalid")
