# -*- coding: utf-8 -*-

from __future__ import absolute_import

from cobra.core.object import Object


class KeyValuePair(Object):
    """
    Class representation of key-value pairs supported in fba-v3

    Parameters
    ----------
    data : dict
        A dictionary containong data about key-value pairs
        {
            "id" : "abc",
            "name" : "abc",
            "key" : "abc",
            "value" : "abc",
            "uri" : "abc"
        }

    Attributes
    ----------
    id : string
        The identifier to associate with this key-value pair
    name : string
        A human readable name for this key-value pair
    key : string
        The key by which we can refer to the value. Must be
        unique in a given list of key-value pair
    value : string
        The value corresponding to that key
    uri : string
        The uri identifies a resource that defines the associated
        key component

    """

    def __init__(self, data={}):
        if "id" in data:
            id = data["id"]
        else:
            id = None
        if "name" in data:
            name = data["name"]
        else:
            name = None
        Object.__init__(self, id, name)
        if "key" in data:
            self._key = data["key"]
        else:
            self._key = None
        if "value" in data:
            self._value = data["value"]
        else:
            self._value = None
        if "uri" in data:
            self._uri = data["uri"]
        else:
            self._uri = None

    @property
    def key(self):
        return getattr(self, "_key", None)

    @key.setter
    def key(self, inKey):
        if not isinstance(inKey, str):
            raise TypeError("the key must be of type string")
        else:
            self._key = inKey

    @property
    def value(self):
        return getattr(self, "_key", None)

    @value.setter
    def value(self, inValue):
        if not isinstance(inValue, str):
            raise TypeError("the value must be of type string")
        else:
            self._value = inValue

    @property
    def uri(self):
        return getattr(self, "_uri", None)

    @uri.setter
    def uri(self, inUri):
        if not isinstance(inUri, str):
            raise TypeError("the uri must be of type string")
        else:
            self._uri = inUri
