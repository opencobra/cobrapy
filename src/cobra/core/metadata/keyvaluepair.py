import collections
from typing import Dict, List, Union


class ListOfKeyValue(collections.MutableSequence):
    """A list extension to store key-value pairs

    Parameters
    ----------
    keyvaluelist : list
        a list containing KeyValueDict objects.
    """

    def __init__(self, keyvaluelist: List = None):
        if keyvaluelist is None:
            keyvaluelist = []
        self._sequence = list()
        if not isinstance(keyvaluelist, list):
            raise TypeError(f"The key-value data must be of type list: {keyvaluelist}")
        else:
            for item in keyvaluelist:
                self.append(item)

    @staticmethod
    def parse_listofkeyvalue(data: List) -> "ListOfKeyValue":
        """Makes an object of ListOfKeyValue using given data"""
        if data is None or isinstance(data, list):
            return ListOfKeyValue(data)
        else:
            raise TypeError(f"Invalid format passed for ListOfKeyValue: {data}")

    def append(self, value: Union[Dict, "KeyValueDict"]) -> None:
        """ Append a KeyValueDict object in its list. """
        if isinstance(value, KeyValueDict):
            self._sequence.append(value)
        elif isinstance(value, dict):
            self._sequence.append(KeyValueDict(**value))
        else:
            raise TypeError(
                f"The data passed for KeyValuePair " f"has invalid format: {value}"
            )

    def insert(self, index: int, value: Union[Dict, "KeyValueDict"]) -> None:
        """Insert a KeyValueDict object at the given index in its list."""
        if isinstance(value, KeyValueDict):
            self._sequence.insert(index, value)
        elif isinstance(value, dict):
            self._sequence.insert(index, KeyValueDict(**value))
        else:
            raise TypeError(
                f"The data passed for KeyValuePair " f"has invalid format: {value}"
            )

    def __getitem__(self, index: int) -> "KeyValueDict":
        return self._sequence[index]

    def __setitem__(self, index: int, value: Union[Dict, "KeyValueDict"]) -> None:
        if isinstance(value, KeyValueDict):
            self._sequence[index] = value
        elif isinstance(value, dict):
            self._sequence[index] = KeyValueDict(**value)
        else:
            raise TypeError(
                f"The data passed for KeyValuePair " f"has invalid format: {value}"
            )

    def __len__(self) -> int:
        return len(self._sequence)

    def __delitem__(self, index) -> None:
        del self._sequence[index]

    def __str__(self) -> str:
        return str(self._sequence)

    def __repr__(self) -> str:
        return f"{self._sequence}"


class KeyValueDict:
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

    def __init__(
        self,
        id: str = None,
        name: str = None,
        key: str = None,
        value: str = None,
        uri: str = None,
    ):
        self.id = id
        self.name = name
        self.key = key
        self.value = value
        self.uri = uri

    @staticmethod
    def parse_keyvaluedict(data: Union[Dict, "KeyValueDict"]) -> "KeyValueDict":
        """ Makes a KeyValueDict object using the data passed. """
        if data is None:
            return KeyValueDict(None)
        elif isinstance(data, KeyValueDict):
            return data
        elif isinstance(data, dict):
            return KeyValueDict(**data)
        else:
            raise TypeError(f"Invalid format for KeyValueDict: '{data}'")

    @property
    def id(self) -> str:
        return self._id

    @id.setter
    def id(self, data: str) -> None:
        if not isinstance(data, str):
            raise TypeError(f"Only string type allowed for 'id': {data}")
        else:
            self._id = data

    @property
    def name(self) -> str:
        return self._name

    @name.setter
    def name(self, data: str) -> None:
        if not isinstance(data, str):
            raise TypeError(f"Only string type allowed for 'name': {data}")
        else:
            self._name = data

    @property
    def key(self) -> str:
        return self._key

    @key.setter
    def key(self, data: str) -> None:
        if not isinstance(data, str):
            raise TypeError(f"Only string type allowed for 'key': {data}")
        else:
            self._key = data

    @property
    def value(self) -> str:
        return self._value

    @value.setter
    def value(self, data: str) -> None:
        if not isinstance(data, str):
            raise TypeError(f"Only string type allowed for 'value': {data}")
        else:
            self._value = data

    @property
    def uri(self) -> str:
        return self._uri

    @uri.setter
    def uri(self, data: str) -> None:
        if not isinstance(data, str):
            raise TypeError(f"Only string type allowed for 'uri': {data}")
        else:
            self._uri = data

    def __str__(self) -> str:
        return str(
            {
                "id": self.id,
                "name": self.name,
                "key": self.key,
                "value": self.value,
                "uri": self.uri,
            }
        )

    def __repr__(self) -> str:
        return self.__str__()
