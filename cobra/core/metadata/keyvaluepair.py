# -*- coding: utf-8 -*-

from __future__ import absolute_import

from collections.abc import MutableMapping, MutableSequence


class ListOfKeyValue(MutableSequence):
    """A list extension to store key-value pairs

    Parameters
    ----------
    creators : list key-value pair data
    """

    def __init__(self, keyvaluelist=None):
        if keyvaluelist is None:
            keyvaluelist = []
        self._sequence = list()
        if not isinstance(keyvaluelist, list):
            raise TypeError("The data passed for ListOfKeyValue "
                            "must be inside a list")
        else:
            for item in keyvaluelist:
                if isinstance(item, self.KeyValuePair):
                    self._sequence.append(item)
                elif isinstance(item, dict):
                    self._sequence.append(self.KeyValuePair(item))
                else:
                    raise TypeError("The data passed for KeyValuepair "
                                    "indexed %s has invalid format"
                                    % keyvaluelist.index(item, 0,
                                                         len(keyvaluelist)))

    def __len__(self):
        return len(self._sequence)

    def __delitem__(self, index):
        del self._sequence[index]

    def insert(self, index, value):
        if isinstance(value, self.KeyValuePair):
            self._sequence.insert(index, value)
        elif isinstance(value, dict):
            self._sequence.insert(index, self.KeyValuePair(value))
        else:
            raise TypeError("The data passed for KeyValuePair "
                            "has invalid format")

    def append(self, value):
        if isinstance(value, self.KeyValuePair):
            self._sequence.append(value)
        elif isinstance(value, dict):
            self._sequence.append(self.KeyValuePair(value))
        else:
            raise TypeError("The data passed for KeyValuePair "
                            "has invalid format")

    def __setitem__(self, index, value):
        if isinstance(value, self.KeyValuePair):
            self._sequence[index] = value
        elif isinstance(value, dict):
            self._sequence[index] = self.KeyValuePair(value)
        else:
            raise TypeError("The data passed for KeyValuePair "
                            "has invalid format")

    def __getitem__(self, index):
        return self._sequence[index]

    def __str__(self):
        return str(self._sequence)

    def __repr__(self):
        return '{}'.format(self._sequence)

    class KeyValuePair(MutableMapping):
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

        def __init__(self, data=None):
            if data is None:
                data = {}
            self._mapping = dict()
            for key, value in data.items():
                if key not in self.VALID_KEYS:
                    raise ValueError("'%s' is not allowed. Only possible "
                                     "keys are : 'id', 'name', 'key', "
                                     "'value', 'uri'" % key)
                if not isinstance(key, str):
                    raise TypeError("All keys must be of type string")
                if not isinstance(value, str):
                    raise TypeError("All values must be of type string")
                self._mapping[key] = value
            for key in self.VALID_KEYS:
                if key not in data:
                    self._mapping[key] = None

        def __getitem__(self, key):
            if key not in self.VALID_KEYS:
                raise ValueError("Key %s is not allowed. Only allowed "
                                 "keys are : 'id', 'name', 'key',"
                                 " 'value', 'uri'" % key)
            return self._mapping[key]

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
