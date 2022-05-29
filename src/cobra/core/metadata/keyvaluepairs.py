from collections import UserDict
from typing import Dict, Iterable, Union
from ...util import format_long_string


class KeyValueEntry:
    """Single key-value entry.

    The key is an attribute on the entry.

    Parameters
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
    def from_data(data: Union[Dict, "KeyValueEntry"]) -> "KeyValueEntry":
        """Makes a KeyValueDict object using the data passed."""
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


class KeyValuePairs(UserDict):
    """An UserDict to store KeyValueEntries.

    Parameters
    ----------
    entries : Iterable
        an iterable containing entry information
    """

    def __init__(self, entries: Iterable[Union[Dict, KeyValueEntry]] = None):
        super().__init__()
        if entries:
            for item in entries:
                entry = KeyValueEntry.from_data(item)
                self.data[entry.key] = entry

    def __setitem__(self, key: str, item: Union[Dict, KeyValueEntry]) -> None:
        entry = KeyValueEntry.from_data(item)
        self.data[key] = entry

    def __str__(self) -> str:
        """Convert KeyValuePairs to str.

        Parameters
        ----------
        self : KeyValuePairs
            UserDict defining key value pairs

        Returns
        ------
        string
            a string representation of a dictionary
        """
        return str(self.to_dict())

    def __repr__(self) -> str:
        return (
            f"{self.__class__.__module__}.{self.__class__.__qualname__}"
            f"({self.to_dict()!r})"
        )

    def _repr_html_(self) -> str:
        return f"""<p><strong>KeyValuePairs</strong></p><p>{format_long_string(
            self.__str__(), 100)}</p>"""

    def to_dict(self) -> dict:
        return {k: v.to_dict() for k, v in self.data.items()}
