from collections import OrderedDict
from typing import Dict, Union, Iterable
from collections.abc import MutableMapping


class KeyValueEntry:
    """ Single key-value entry.

    The key is an attribute on the entry.

    Parameters
    ----------
        id : str
        name : str
        key : str
        value : str
        uri : str
    """
    def __init__(self, id: str = None, name: str = None,
                 key: str = None, value: str = None, uri: str = None):
        self.id = id
        self.name = name
        self.key = key
        self.value = value
        self.uri = uri

    @staticmethod
    def from_data(data: Union[Dict, 'KeyValueEntry']) -> 'KeyValueEntry':
        """ Makes a KeyValueDict object using the data passed. """
        if data is None:
            return KeyValueEntry()
        elif isinstance(data, KeyValueEntry):
            return data
        elif isinstance(data, dict):
            return KeyValueEntry(**data)
        else:
            raise TypeError(f"Invalid format for KeyValueEntry: '{data}'")

    def to_dict(self) -> dict:
        return {
            "id": self.id,
            "name": self.name,
            "key": self.key,
            "value": self.value,
            "uri": self.uri,
        }

    def __str__(self) -> str:
        return str(self.to_dict())

    def __repr__(self) -> str:
        return f"<{self.__class__.__name__} ({self.key}, {self.value}, {self.uri})>"


class KeyValuePairs(MutableMapping):
    """An OrderedDict extension to store KeyValueEntries

    Parameters
    ----------
    entries : Iterable
        an iterable containing entry information
    """
    def __init__(self, entries: Iterable[Union[Dict, KeyValueEntry]] = None):
        self.mapping = OrderedDict()  # type: OrderedDict[str, KeyValueEntry]
        if entries:
            for item in entries:
                entry = KeyValueEntry.from_data(item)
                self.mapping[entry.key] = entry

    def __getitem__(self, key: str) -> KeyValueEntry:
        return self.mapping.__getitem__(key)

    def __setitem__(self, key: str, item: Union[Dict, KeyValueEntry]) -> None:
        entry = KeyValueEntry.from_data(item)
        self.mapping[key] = entry

    def __len__(self) -> int:
        return len(self.mapping)

    def __iter__(self):
        return iter(self.mapping)

    def __delitem__(self, key: str) -> None:
        del self.mapping[key]

    def __str__(self) -> str:
        return str(self.mapping)

    def __repr__(self) -> str:
        return f"<{self.__class__.__name__} [{len(self)}]>"
