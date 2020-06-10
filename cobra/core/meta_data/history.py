# -*- coding: utf-8 -*-

from __future__ import absolute_import

import datetime


# The possible keys inside creator dict
CREATOR_KEYS = ["first_name", "last_name", "email", "organization_name"]


def validateDate(date_text):
    """Validate if the date format is of type w3cdtf ISO 8601"""
    try:
        datetime.datetime.strptime(date_text, '%Y-%m-%dT%H:%M:%S%z')
    except ValueError as e:
        raise ValueError(str(e))
    return True


class History:
    """
    Class representation of history of a given component i.e. creator,
    created date and modification dates

    Parameters
    ----------
    creator : dict
        A dictionary containong details of creator's name, email and
        organisation name
    created : string
        The date when component is created in W3CDTF ISO 8601 format
    modified : list
        A list of dates about the component modification

    Attributes
    ----------
    creator : dict
        A dictionary containong details of creator's name, email and
        organisation name
    created : string
        The date when component is created in W3CDTF ISO 8601 format
    modified : list
        A list of dates about the component modification

    """

    def __init__(self, creators=[], created=None, modified=[]):
        if creators is None:
            creators = []
        self._creators = self.ListOfCreators(creators)
        if isinstance(created, str):
            validateDate(created)
            self._created = created
        elif created is None:
            self._created = None
        else:
            raise TypeError('Only None and string types are possible for '
                            '"created" date attribute')
        if modified is None:
            modified = []
        self._modified = self.ModifiedHistory(modified)

    @property
    def creators(self):
        return self._creators

    @creators.setter
    def creators(self, creators_list):
        self._creators = self.ListOfCreators(creators_list)

    @property
    def created(self):
        return self._created

    @created.setter
    def created(self, value):
        if not isinstance(value, str):
            raise TypeError("The date passed must be a string")
        else:
            validateDate(value)
            self._created = value

    @property
    def modified(self):
        return self._modified

    @modified.setter
    def modified(self, value):
        self._modified = self.ModifiedHistory(value)

    class ListOfCreators(list):
        """A list extension to store each creator's info

        Parameters
        ----------
        creators : list containing info about creators
        """

        def __init__(self, creators=[]):
            if not isinstance(creators, list):
                raise TypeError("The data passed for creators must be "
                                "inside a list")
            else:
                for item in creators:
                    if isinstance(item, History.Creator):
                        list.append(self, item)
                    elif isinstance(item, dict):
                        list.append(self, History.Creator(item))
                    else:
                        raise TypeError("The data passed for creator "
                                        "indexed %s has invalid format"
                                        % creators.index(item, 0,
                                                         len(creators)))

        def __len__(self):
            return list.__len__(self)

        def __delitem__(self, index):
            list.__delitem__(self, index)

        def insert(self, index, value):
            if isinstance(value, History.Creator):
                list.insert(self, index, value)
            elif isinstance(value, dict):
                list.insert(self, index, History.Creator(value))
            else:
                raise TypeError("The data passed has invalid format")

        def append(self, value):
            if isinstance(value, History.Creator):
                list.append(self, value)
            elif isinstance(value, dict):
                list.append(self, History.Creator(value))
            else:
                raise TypeError("The data passed has invalid format")

        def __setitem__(self, index, value):
            if isinstance(value, History.Creator):
                list.__setitem__(self, index, value)
            elif isinstance(value, dict):
                list.__setitem__(self, index, History.Creator(value))
            else:
                raise TypeError("The data passed has invalid format")

        def __getitem__(self, index):
            return list.__getitem__(self, index)

    class Creator(dict):
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

        def __init__(self, creator_dict={}):
            if not isinstance(creator_dict, dict):
                raise TypeError("The value passed for creator must "
                                "be of type dict.")
            for key in CREATOR_KEYS:
                if key not in creator_dict:
                    dict.__setitem__(self, key, None)
                else:
                    if not isinstance(creator_dict[key], str):
                        raise TypeError("All the values passed must "
                                        "be of type string.")
                    else:
                        dict.__setitem__(self, key, creator_dict[key])

        def __getitem__(self, key):
            if key not in CREATOR_KEYS:
                raise ValueError("Key %s is not allowed. only allowed "
                                 "keys are 'first_name', 'last_name', "
                                 "'email', 'organization_name'." % key)
            return dict.__getitem__(self, key)

        def __setitem__(self, key, value):
            if key not in CREATOR_KEYS:
                raise ValueError("Key %s is not allowed. Only allowed "
                                 "keys are 'first_name', 'last_name', "
                                 "'email', 'organization_name'." % key)
            if not isinstance(value, str):
                raise TypeError("Value passed must be of type string.")
            dict.__setitem__(self, key, value)

        def __delitem__(self, key):
            dict.__delitem__(self, key)

        def __iter__(self):
            return dict.__iter__(self)

        def __len__(self):
            return dict.__len__(self)

        def __contains__(self, x):
            return dict.__contains__(self, x)

    class ModifiedHistory(list):
        """A list extension to store modification dates. Only Restricted
        type of entries are possible.

        Parameters
        ----------
        modifiedList : list containing modification dates in W3CDTF ISO
                       8601 format
        """

        def __init__(self, modifiedList=[]):
            if not isinstance(modifiedList, list):
                raise TypeError("The dates passed must be inside a list")
            for item in modifiedList:
                if not isinstance(item, str):
                    raise ValueError("Each date must be of type string")
                validateDate(item)
                list.append(self, item)

        def __len__(self):
            return list.__len__(self)

        def __delitem__(self, index):
            list.__delitem__(self, index)

        def insert(self, index, value):
            validateDate(value)
            list.insert(self, index, value)

        def append(self, value):
            validateDate(value)
            list.append(self, value)

        def __setitem__(self, index, value):
            validateDate(value)
            list.__setitem__(self, index, value)

        def __getitem__(self, index):
            return list.__getitem__(self, index)
