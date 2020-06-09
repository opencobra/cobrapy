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

    def __init__(self, creator={}, created=None, modified=[]):
        self._creator = self.Creator(creator)
        if isinstance(created, str):
            validateDate(created)
            self._created = created
        elif created == None:
            self._created = None
        else:
            raise TypeError("Only None and string type possible for created date")
        self._modified = self.ModifiedHistory(modified)

    @property
    def creator(self):
        return self._creator

    @creator.setter
    def creator(self, creator_dict):
        self._creator = self.Creator(creator_dict)

    @property
    def created(self):
        return self._created

    @created.setter
    def created(self, value):
        if not isinstance(value, str):
            raise TypeError("date passed must be a string")
        else:
            validateDate(value)
            self._created = value

    @property
    def modified(self):
        return self._modified

    @modified.setter
    def modified(self, value):
        self._modified = self.ModifiedHistory(value)


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
                raise TypeError("The value passed must be of type dict.")
            for key in CREATOR_KEYS:
                if key not in creator_dict:
                    dict.__setitem__(self, key, None)
                else:
                    if not isinstance(creator_dict[key], str):
                        raise TypeError("All the values passed must be of type string.")
                    else:
                        dict.__setitem__(self, key, creator_dict[key])

        def __getitem__(self,key):
            if key not in CREATOR_KEYS:
                raise ValueError("Key %s is not allowed. Allowed keys are 'first_name', 'last_name', 'email', 'organization_name'." %key)
            return dict.__getitem__(self,key)

        def __setitem__(self, key, value):
            if key not in CREATOR_KEYS:
                raise ValueError("Key %s is not allowed. Allowed keys are 'first_name', 'last_name', 'email', 'organization_name'." %key)
            if not isinstance(value, str):
                raise TypeError("Value passed must be of type string.")
            dict.__setitem__(self,key,value)

        def __delitem__(self, key):
            dict.__delitem__(self,key)

        def __iter__(self):
            return dict.__iter__(self)

        def __len__(self):
            return dict.__len__(self)

        def __contains__(self, x):
            return dict.__contains__(self,x)


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
            list.__delitem__(self,index)

        def insert(self, index, value):
            validateDate(value)
            list.insert(self,index, value)

        def append(self, value):
            validateDate(value)
            list.append(self,value)

        def __setitem__(self, index, value):
            validateDate(value)
            list.__setitem__(self,index, value)

        def __getitem__(self, index):
            return list.__getitem__(self,index)
