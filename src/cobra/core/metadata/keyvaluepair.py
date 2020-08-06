# -*- coding: utf-8 -*-

import collections

try:
    collectionsAbc = collections.abc
except AttributeError:
    collectionsAbc = collections


class ListOfKeyValue(collectionsAbc.MutableSequence):
    """A list extension to store key-value pairs

    Parameters
    ----------
    keyvaluelist : list
        a list containing KeyValueDict objects.
    """

    def __init__(self, keyvaluelist: 'list' = None):
        if keyvaluelist is None:
            keyvaluelist = []
        self._sequence = list()
        if not isinstance(keyvaluelist, list):
            raise TypeError("The key-value data must be of"
                            " type list: {}".format(keyvaluelist))
        else:
            for item in keyvaluelist:
                self.append(item)

    def parse_listofKeyValue(data):
        if data is None or isinstance(data, list):
            return ListOfKeyValue(data)
        else:
            raise TypeError("Invalid format passed "
                            "for ListOfKeyValue: {}".format(data))

    def append(self, value):
        """ Append a KeyValueDict object in its list.
        """
        if isinstance(value, KeyValueDict):
            self._sequence.append(value)
        elif isinstance(value, dict):
            self._sequence.append(KeyValueDict(**value))
        else:
            raise TypeError("The data passed for KeyValuePair "
                            "has invalid format: {}".format(value))

    def insert(self, index, value):
        """ Insert a KeyValueDict object at
        the given index in its list.
        """
        if isinstance(value, KeyValueDict):
            self._sequence.insert(index, value)
        elif isinstance(value, dict):
            self._sequence.insert(index, KeyValueDict(**value))
        else:
            raise TypeError("The data passed for KeyValuePair "
                            "has invalid format: {}".format(value))

    def __getitem__(self, index):
        return self._sequence[index]

    def __setitem__(self, index, value):
        if isinstance(value, KeyValueDict):
            self._sequence[index] = value
        elif isinstance(value, dict):
            self._sequence[index] = KeyValueDict(**value)
        else:
            raise TypeError("The data passed for KeyValuePair "
                            "has invalid format: {}".format(value))

    def __len__(self):
        return len(self._sequence)

    def __delitem__(self, index):
        del self._sequence[index]

    def __str__(self):
        return str(self._sequence)

    def __repr__(self):
        return '{}'.format(self._sequence)


class KeyValueDict(object):
    """ Class representation of a single key-value data
    for fbc-v3 key-value pair data.
    Parameters
    ----------
        id : str
        name : str
        key : str
        value : str
        uri : str
    Attributes
    ----------
        id : str
        name : str
        key : str
        value : str
        uri : str
    """

    def __init__(self, id: str = None, name: str = None, key: str = None,
                 value: str = None, uri: str = None):
        self._id = None
        self._name = None
        self._key = None
        self._value = None
        self._uri = None
        self.id = id
        self.name = name
        self.key = key
        self.value = value
        self.uri = uri

    @staticmethod
    def parse_keyValueDict(data) -> 'KeyValueDict':
        """Tries to parse the KeyValueDict."""
        if data is None:
            return KeyValueDict(None)
        elif isinstance(data, KeyValueDict):
            return data
        elif isinstance(data, dict):
            return KeyValueDict(**data)
        else:
            raise TypeError("Invalid format for KeyValueDict:"
                            " '{}'".format(data))

    @property
    def id(self):
        return getattr(self, "_id", None)

    @id.setter
    def id(self, data):
        if not isinstance(data, str):
            raise TypeError("Only string type allowed "
                            "for 'id': {}".format(data))
        else:
            self._id = data

    @property
    def name(self):
        return getattr(self, "_name", None)

    @name.setter
    def name(self, data):
        if not isinstance(data, str):
            raise TypeError("Only string type allowed "
                            "for 'name': {}".format(data))
        else:
            self._name = data

    @property
    def key(self):
        return getattr(self, "_key", None)

    @key.setter
    def key(self, data):
        if not isinstance(data, str):
            raise TypeError("Only string type allowed "
                            "for 'key': {}".format(data))
        else:
            self._key = data

    @property
    def value(self):
        return getattr(self, "_value", None)

    @value.setter
    def value(self, data):
        if not isinstance(data, str):
            raise TypeError("Only string type allowed "
                            "for 'value': {}".format(data))
        else:
            self._value = data

    @property
    def uri(self):
        return getattr(self, "_uri", None)

    @uri.setter
    def uri(self, data):
        if not isinstance(data, str):
            raise TypeError("Only string type allowed for "
                            "'uri': {}".format(data))
        else:
            self._uri = data

    def __str__(self):
        return str({
            "id": self.id,
            "name": self.name,
            "key": self.key,
            "value": self.value,
            "uri": self.uri})

    def __repr__(self):
        return self.__str__()
