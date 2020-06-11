# -*- coding: utf-8 -*-

from __future__ import absolute_import


class ListOfKeyValue(list):
    """A list extension to store key-value pairs

    Parameters
    ----------
    creators : list key-value pair data
    """

    def __init__(self, keyvaluelist=[]):
        if not isinstance(keyvaluelist, list):
            raise TypeError("The data passed for ListOfKeyValue "
                            "must be inside a list")
        else:
            for item in keyvaluelist:
                if isinstance(item, self.KeyValuePair):
                    list.append(self, item)
                elif isinstance(item, dict):
                    list.append(self, self.KeyValuePair(item))
                else:
                    raise TypeError("The data passed for KeyValuepair "
                                    "indexed %s has invalid format"
                                    % keyvaluelist.index(item, 0,
                                                         len(keyvaluelist)))

    def __len__(self):
        return list.__len__(self)

    def __delitem__(self, index):
        list.__delitem__(self, index)

    def insert(self, index, value):
        if isinstance(value, self.KeyValuePair):
            list.insert(self, index, value)
        elif isinstance(value, dict):
            list.insert(self, index, self.KeyValuePair(value))
        else:
            raise TypeError("The data passed for KeyValuePair "
                            "has invalid format")

    def append(self, value):
        if isinstance(value, self.KeyValuePair):
            list.append(self, value)
        elif isinstance(value, dict):
            list.append(self, self.KeyValuePair(value))
        else:
            raise TypeError("The data passed for KeyValuePair "
                            "has invalid format")

    def __setitem__(self, index, value):
        if isinstance(value, self.KeyValuePair):
            list.__setitem__(self, index, value)
        elif isinstance(value, dict):
            list.__setitem__(self, index, self.KeyValuePair(value))
        else:
            raise TypeError("The data passed for KeyValuePair "
                            "has invalid format")

    def __getitem__(self, index):
        return list.__getitem__(self, index)

    class KeyValuePair(dict):
        """
        Class representation of key-value pairs supported in fba-v3

        Parameters
        ----------
        data : dict
            A dictionary containing data about key-value pairs
            {
                "id" : "abc",
                "name" : "abc",
                "key" : "abc",
                "value" : "abc",
                "uri" : "abc"
            }

        """

        VALID_KEYS = ["id", "name", "key", "value", "uri"]

        def __init__(self, data={}):
            for key, value in data.items():
                if key not in self.VALID_KEYS:
                    raise ValueError("'%s' is not allowed. Only possible "
                                     "keys are : 'id', 'name', 'key', "
                                     "'value', 'uri'" % key)
                if not isinstance(key, str):
                    raise TypeError("All keys must be of type string")
                if not isinstance(value, str):
                    raise TypeError("All values must be of type string")
                dict.__setitem__(self, key, value)
            for key in self.VALID_KEYS:
                if key not in data:
                    dict.__setitem__(self, key, None)

        def __getitem__(self, key):
            if key not in self.VALID_KEYS:
                raise ValueError("Key %s is not allowed. Only allowed "
                                 "keys are : 'id', 'name', 'key',"
                                 " 'value', 'uri'" % key)
            return dict.__getitem__(self, key)

        def __setitem__(self, key, value):
            """Restricting the keys and values that can be set.
               Only allowed keys are : 'id', 'name', 'key', 'value', 'uri''
            """
            if key not in self.VALID_KEYS:
                raise ValueError("Key %s is not allowed. Only allowed"
                                 " keys are : 'id', 'name', 'key',"
                                 " 'value', 'uri'" % key)
            if not isinstance(value, str):
                raise TypeError("The value must be of type string")
            dict.__setitem__(self, key, value)

        def __delitem__(self, key):
            dict.__delitem__(self, key)

        def __iter__(self):
            return dict.__iter__(self)

        def __len__(self):
            return dict.__len__(self)

        def __contains__(self, x):
            return dict.__contains__(self, x)
