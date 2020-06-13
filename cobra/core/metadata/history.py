# -*- coding: utf-8 -*-

from __future__ import absolute_import

import datetime
from collections.abc import MutableMapping, MutableSequence


# The possible keys inside creator dict
CREATOR_KEYS = ["first_name", "last_name", "email", "organization_name"]


def validateDate(date_text):
    """Validate if the date format is of type w3cdtf ISO 8601"""
    try:
        datetime.datetime.strptime(date_text, '%Y-%m-%dT%H:%M:%S%z')
    except ValueError as e:
        raise ValueError(str(e))
    return True


class History(MutableMapping):
    """
    Class representation of history of a given component i.e. creator,
    created date and modification dates. It is basically an extended
    dictionary with some restrictions

    Parameters
    ----------
    creator : dict
        A dictionary containong details of creator's name, email and
        organisation name
    created : string
        The date when component is created in W3CDTF ISO 8601 format
    modified : list
        A list of dates about the component modification

    Allowed Keys
    ----------
    creator : dict
        A dictionary containong details of creator's name, email and
        organisation name
    created : string
        The date when component is created in W3CDTF ISO 8601 format
    modified : list
        A list of dates about the component modification

    """

    VALID_KEYS = ["creators", "created", "modified"]

    def __init__(self, creators=None, created=None, modified=None):
        if creators is None:
            creators = []
        if modified is None:
            modified = []
        self._mapping = dict()
        self._mapping["creators"] = self.ListOfCreators(creators)
        if isinstance(created, str):
            validateDate(created)
            self._mapping["created"] = created
        elif created is None:
            self._mapping["created"] = None
        else:
            raise TypeError('Only None and string types are possible for '
                            '"created" date attribute')
        self._mapping["modified"] = self.ModifiedHistory(modified)

    def __getitem__(self, key):
        if key not in self.VALID_KEYS:
            raise ValueError("Key %s is not allowed. Only allowed "
                             "keys are : 'creators', 'created', 'modified'"
                             % key)
        return self._mapping[key]

    def __setitem__(self, key, value):
        """Restricting the keys and values that can be set.
           Only allowed keys are : 'id', 'name', 'key', 'value', 'uri''
        """
        if key not in self.VALID_KEYS:
            raise ValueError("Key %s is not allowed. Only allowed"
                             " keys are : 'id', 'name', 'key',"
                             " 'value', 'uri'" % key)
        if key == "creators":
            if isinstance(value, self.ListOfCreators):
                self._mapping[key] = value
            elif isinstance(value, list):
                self._mapping[key] = self.ListOfCreators(value)
            else:
                raise TypeError("The passed format for creators is invalid")
        elif key == "created":
            if not isinstance(value, str):
                raise TypeError("The date passed must be a string")
            else:
                validateDate(value)
                self._mapping[key] = value
        elif key == "modified":
            if isinstance(value, self.ModifiedHistory):
                self._mapping[key] = value
            elif isinstance(value, list):
                self._mapping[key] = self.ModifiedHistory(value)
            else:
                raise TypeError("The passed format for modification"
                                " history is invalid")

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

    class ListOfCreators(MutableSequence):
        """A list extension to store each creator's info

        Parameters
        ----------
        creators : list containing info about creators
        """

        def __init__(self, creators=None):
            if creators is None:
                creators = []
            self._sequence = list()
            if not isinstance(creators, list):
                raise TypeError("The data passed for creators must be "
                                "inside a list")
            else:
                for item in creators:
                    if isinstance(item, History.Creator):
                        self._sequence.append(item)
                    elif isinstance(item, dict):
                        self._sequence.append(History.Creator(item))
                    else:
                        raise TypeError("The data passed for creator "
                                        "indexed %s has invalid format"
                                        % creators.index(item, 0,
                                                         len(creators)))

        def __len__(self):
            return len(self._sequence)

        def __delitem__(self, index):
            del self._sequence[index]

        def insert(self, index, value):
            if isinstance(value, History.Creator):
                self._sequence.insert(index, value)
            elif isinstance(value, dict):
                self._sequence.insert(index, History.Creator(value))
            else:
                raise TypeError("The data passed has invalid format")

        def append(self, value):
            if isinstance(value, History.Creator):
                self._sequence.append(value)
            elif isinstance(value, dict):
                self._sequence.append(History.Creator(value))
            else:
                raise TypeError("The data passed has invalid format")

        def __setitem__(self, index, value):
            if isinstance(value, History.Creator):
                self._sequence[index] = value
            elif isinstance(value, dict):
                self._sequence[index] = History.Creator(value)
            else:
                raise TypeError("The data passed has invalid format")

        def __getitem__(self, index):
            return self._sequence[index]

        def __str__(self):
            return str(self._sequence)

        def __repr__(self):
            return '{}'.format(self._sequence)

    class Creator(MutableMapping):
        """A dictionary extension to store basic info of this component
           creator

        Parameters
        ----------
        creator_dict : dict containing info about creator
            {
                "first_name" : "abc",
                "last_name" : "abc",
                "email" : "abc",
                "organization_name" : "abc"
            }
        """

        def __init__(self, creator_dict=None):
            if creator_dict is None:
                creator_dict = {}
            self._mapping = dict()
            if not isinstance(creator_dict, dict):
                raise TypeError("The value passed for creator must "
                                "be of type dict.")
            for key in CREATOR_KEYS:
                if key not in creator_dict:
                    self._mapping[key] = None
                else:
                    if not isinstance(creator_dict[key], str):
                        raise TypeError("All the values passed must "
                                        "be of type string.")
                    else:
                        self._mapping[key] = creator_dict[key]

        def __getitem__(self, key):
            if key not in CREATOR_KEYS:
                raise ValueError("Key %s is not allowed. only allowed "
                                 "keys are 'first_name', 'last_name', "
                                 "'email', 'organization_name'." % key)
            return self._mapping[key]

        def __setitem__(self, key, value):
            if key not in CREATOR_KEYS:
                raise ValueError("Key %s is not allowed. Only allowed "
                                 "keys are 'first_name', 'last_name', "
                                 "'email', 'organization_name'." % key)
            if not isinstance(value, str):
                raise TypeError("Value passed must be of type string.")
            self._mapping[key] = value

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

    class ModifiedHistory(MutableSequence):
        """A list extension to store modification dates. Only Restricted
        type of entries are possible.

        Parameters
        ----------
        modifiedList : list containing modification dates in W3CDTF ISO
                       8601 format
        """

        def __init__(self, modifiedList=None):
            if modifiedList is None:
                modifiedList = []
            self._sequence = list()
            if not isinstance(modifiedList, list):
                raise TypeError("The dates passed must be inside a list")
            for item in modifiedList:
                if not isinstance(item, str):
                    raise ValueError("Each date must be of type string")
                validateDate(item)
                self._sequence.append(item)

        def __len__(self):
            return len(self._sequence)

        def __delitem__(self, index):
            del self._sequence[index]

        def insert(self, index, value):
            validateDate(value)
            self._sequence.insert(index, value)

        def append(self, value):
            validateDate(value)
            self._sequence.append(value)

        def __setitem__(self, index, value):
            validateDate(value)
            self._sequence[index] = value

        def __getitem__(self, index):
            return self._sequence[index]

        def __str__(self):
            return str(self._sequence)

        def __repr__(self):
            return '{}'.format(self._sequence)
