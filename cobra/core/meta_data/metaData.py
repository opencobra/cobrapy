# -*- coding: utf-8 -*-

from __future__ import absolute_import

from cobra.core.meta_data.cvTerm import CVTerm
from cobra.core.meta_data.history import History
from cobra.core.meta_data.keyValuePair import ListOfKeyValue


class MetaData(dict):
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

    VALID_KEYS = ["sbo", "cvTerms", "history", "listofKeyValue"]

    def __init__(self, cvterm=None, history=None, listofKeyValue=None):
        # setting the cvterm
        if cvterm is None:
            dict.__setitem__(self, "cvTerms", CVTerm())
        elif isinstance(cvterm, CVTerm):
            dict.__setitem__(self, "cvTerms", cvterm)
        elif isinstance(cvterm, dict):
            dict.__setitem__(self, "cvTerms", CVTerm(cvterm))
        else:
            raise TypeError("Invalid format passed for cvterm")
        # setting the history of the component
        if history is None:
            dict.__setitem__(self, "history", History())
        elif isinstance(history, History):
            dict.__setitem__(self, "history", history)
        elif isinstance(history, dict):
            if "creators" not in history:
                history["creators"] = []
            if "created" not in history:
                history["created"] = None
            if "modified" not in history:
                history["modified"] = []
            dict.__setitem__(self, "history", History(history["creators"],
                                                      history["created"],
                                                      history["modified"]))
        else:
            raise TypeError("Invalid format passed for history")
        # setting the list of key-value pair
        if listofKeyValue is not None:
            if isinstance(listofKeyValue, ListOfKeyValue):
                dict.__setitem__(self, "listofKeyValue", listofKeyValue)
            elif isinstance(listofKeyValue, list):
                dict.__setitem__(self, "listofKeyValue",
                                 ListOfKeyValue(listofKeyValue))
            else:
                raise TypeError("Key value pairs must be passed in a list")
        else:
            dict.__setitem__(self, "listofKeyValue", ListOfKeyValue())

    def __getitem__(self, key):
        if key not in self.VALID_KEYS:
            raise ValueError("Key %s is not allowed. Only allowed "
                             "keys are : 'sbo', 'cvTerms', 'history', "
                             "'listofKeyValue'" % key)
        return dict.__getitem__(self, key)

    def __setitem__(self, key, value):
        """Restricting the keys and values that can be set.
           Only allowed keys are : 'sbo', 'cvTerms', 'history',
           'listofKeyValue'
        """
        if key not in self.VALID_KEYS:
            raise ValueError("Key %s is not allowed. Only allowed "
                             "keys are : 'sbo', 'cvTerms', 'history', "
                             "'listofKeyValue'" % key)
        if key == "cvTerms":
            if isinstance(value, CVTerm):
                dict.__setitem__(self, "cvTerms", value)
            elif isinstance(value, dict):
                dict.__setitem__(self, "cvTerms", CVTerm(value))
            else:
                raise TypeError("This passed format for cvterm is invalid")
        elif key == "history":
            if isinstance(history, History):
                dict.__setitem__(self, "history", history)
            elif isinstance(history, dict):
                if "creators" not in history:
                    history["creators"] = []
                if "created" not in history:
                    history["created"] = None
                if "modified" not in history:
                    history["modified"] = []
                dict.__setitem__(self, "history", History(history["creators"],
                                                          history["created"],
                                                          history["modified"]))
        elif key == "listofKeyValue":
            if isinstance(listofKeyValue, ListOfKeyValue):
                dict.__setitem__(self, "listofKeyValue", listofKeyValue)
            elif isinstance(listofKeyValue, list):
                dict.__setitem__(self, "listofKeyValue",
                                 ListOfKeyValue(listofKeyValue))
            else:
                raise TypeError("Key value pairs must be passed in a list")
        elif key == "sbo":
            dict.__setitem__(self, "sbo", value)

    def __delitem__(self, key):
        dict.__delitem__(self, key)

    def __iter__(self):
        return dict.__iter__(self)

    def __len__(self):
        return dict.__len__(self)

    def __contains__(self, x):
        return dict.__contains__(self, x)
