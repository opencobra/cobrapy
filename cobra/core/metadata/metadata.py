# -*- coding: utf-8 -*-

from __future__ import absolute_import

from collections import defaultdict
from collections.abc import MutableMapping

from cobra.core.metadata.cvterm import CVTerm
from cobra.core.metadata.history import History
from cobra.core.metadata.keyvaluepair import ListOfKeyValue


class MetaData(MutableMapping):
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
    listofkeyvalue : list
        Some key-value pairs which are not suitable to be
        represented anywhere else in the model.

    """

    VALID_KEYS = ["sbo", "cvterms", "history", "listofkeyvalue", "annotation"]

    def __init__(self, cvterm=None, history=None, listofkeyvalue=None):
        self._mapping = dict()
        self._mapping["annotation"] = defaultdict(list)
        # setting the cvterm
        if cvterm is None:
            self._mapping["cvterms"] = CVTerm(self, None)
        elif isinstance(cvterm, CVTerm):
            self._mapping["cvterms"] = cvterm
        elif isinstance(cvterm, dict):
            self._mapping["cvterms"] = CVTerm(self, cvterm)
        else:
            raise TypeError("Invalid format passed for cvterm")
        # setting the history of the component
        if history is None:
            self._mapping["history"] = History()
        elif isinstance(history, History):
            self._mapping["history"] = history
        elif isinstance(history, dict):
            if "creators" not in history:
                history["creators"] = []
            if "created" not in history:
                history["created"] = None
            if "modified" not in history:
                history["modified"] = []
            self._mapping["history"] = History(history["creators"],
                                               history["created"],
                                               history["modified"])
        else:
            raise TypeError("Invalid format passed for history")
        # setting the list of key-value pair
        if listofkeyvalue is not None:
            if isinstance(listofkeyvalue, ListOfKeyValue):
                self._mapping["listofkeyvalue"] = listofkeyvalue
            elif isinstance(listofkeyvalue, list):
                self._mapping["listofkeyvalue"] = ListOfKeyValue(listofkeyvalue)
            else:
                raise TypeError("Key value pairs must be passed in a list")
        else:
            self._mapping["listofkeyvalue"] = ListOfKeyValue()

    def __getitem__(self, key):
        if key not in self.VALID_KEYS:
            raise ValueError("Key %s is not allowed. Only allowed "
                             "keys are : 'sbo', 'cvterms', 'history', "
                             "'listofkeyvalue', 'annotation'" % key)
        return self._mapping[key]

    def __setitem__(self, key, value):
        """Restricting the keys and values that can be set.
           Only allowed keys are : 'sbo', 'cvterms', 'history',
           'listofkeyvalue'
        """
        if key not in self.VALID_KEYS:
            raise ValueError("Key %s is not allowed. Only allowed "
                             "keys are : 'sbo', 'cvterms', 'history', "
                             "'listofkeyvalue', 'annotation'" % key)
        if key == "cvterms":
            if isinstance(value, CVTerm):
                self._mapping["cvterms"] = value
            elif isinstance(value, dict):
                self._mapping["cvterms"] = CVTerm(self, value)
            else:
                raise TypeError("This passed format for cvterm is invalid")
        elif key == "history":
            if isinstance(history, History):
                self._mapping["history"] = history
            elif isinstance(history, dict):
                if "creators" not in history:
                    history["creators"] = []
                if "created" not in history:
                    history["created"] = None
                if "modified" not in history:
                    history["modified"] = []
                self._mapping["history"] = History(history["creators"],
                                                   history["created"],
                                                   history["modified"])
        elif key == "listofkeyvalue":
            if isinstance(listofkeyvalue, ListOfKeyValue):
                self._mapping["listofkeyvalue"] = listofkeyvalue
            elif isinstance(listofkeyvalue, list):
                self._mapping["listofkeyvalue"] = ListOfKeyValue(listofkeyvalue)
            else:
                raise TypeError("Key value pairs must be passed in a list")
        elif key == "sbo":
            self._mapping["sbo"] = value
        elif key == "annotation":
            self._mapping["annotation"] = value

    def __delitem__(self, key):
        del self._mapping[key]

    def __iter__(self):
        return iter(self._mapping)

    def __len__(self):
        return len(self._mapping)

    def __str__(self):
        return str(self._mapping)

    def __repr__(self):
        return '{}'.format(self._mapping)
