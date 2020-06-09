# -*- coding: utf-8 -*-

from __future__ import absolute_import

from cobra.core.meta_data.cvTerm import CVTerm
from cobra.core.meta_data.history import History
from cobra.core.meta_data.keyValuePair import KeyValuePair


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
        if isinstance(cvterm, CVTerm):
            self._cvTerms = cvterm
        elif isinstance(cvterm, dict):
            self._cvTerms = CVTerm(cvterm)
        else:
            self._cvTerms = CVTerm()
        # setting the history of the component
        if isinstance(history, History):
            self._history = history
        elif isinstance(history, dict):
            if "creator" not in history:
                history["creator"] = {}
            if "created" not in history:
                history["created"] = None
            if "modified" not in history:
                history["modified"] = []
            self._history = History(history["creator"], history["created"], history["modified"])
        else:
            self._history = History()
        # setting the list of key-value pair
        if listofKeyValue is not None:
            if isinstance(listofKeyValue, list):
                keyValue = []
                for item in listofKeyValue:
                    if isinstance(item, dict):
                        keyValue.append(KeyValuePair(item))
                    if isinstance(item, KeyValuePair):
                        keyValue.append(item)
                    else:
                        raise TypeError("Each entry of key-value pair must be a dict or KeyValuePair object")
                self.listofKeyValue = keyValue
            else:
                raise TypeError("Key value pairs must be passed in a list")
        else:
            self.listofKeyValue = []

    @property
    def cvTerms(self):
        return getattr(self, "_cvTerms", None)

    @cvTerms.setter
    def cvTerms(self, value):
        if value == self._cvTerms:
            pass
        elif isinstance(value, CVTerm):
            self._cvTerms = value
        elif isinstance(cvterm, dict):
            self._cvTerms = CVTerm(cvterm)
        else:
            raise TypeError("This passed format for cvterm is not acceptable")

    @property
    def history(self):
        return getattr(self, "_history", None)

    @history.setter
    def history(self, value):
        if value == self._history:
            pass
        elif isinstance(value, History):
            self._history = value
        elif isinstance(history, dict):
            if "creator" not in history:
                history["creator"] = {}
            if "created" not in history:
                history["created"] = None
            if "modified" not in history:
                history["modified"] = []
            self._history = History(value["creator"], value["created"], value["modified"])
        else:
            raise TypeError("This passed format for history is not acceptable")
