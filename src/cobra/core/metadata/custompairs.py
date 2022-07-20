"""Class to deal with Key-Value pairs.

Key-Value pairs are described in SBML FBC3. For the FBC3 standard, see
https://github.com/bgoli/sbml-fbc-spec/blob/main/sf_svn/spec/main.pdf
"""

from collections import UserDict
from typing import Dict, Iterable, Optional, Union

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
        id: Optional[str] = None,
        name: Optional[str] = None,
        key: Optional[str] = None,
        value: Optional[str] = None,
        uri: Optional[str] = None,
    ):
        """Initialize the KeyValueEntry class.

        Parameters
        ----------
        id: str
            optional. Default None.
        name: str
            optional. Default None.
        key: str
            optional. Default None.
        value: str
            optional. Default None.
        uri: str
            optional. Default None.
        """
        self.id = id
        self.name = name
        self.key = key
        self.value = value
        self.uri = uri

    @staticmethod
    def from_data(data: Optional[Union[Dict, "KeyValueEntry"]]) -> "KeyValueEntry":
        """Make a KeyValueDict object using the data passed.

        Parameters
        ----------
        data - dict or KeyValueEntry
            If dict, will use the values of the dictionary to populate a new
            KeyValueEntry. If None, will return empty KeyValueEntry.

        Returns
        -------
        KeyValueEntry
        """
        if data is None:
            return KeyValueEntry()
        elif isinstance(data, KeyValueEntry):
            return data
        elif isinstance(data, dict):
            return KeyValueEntry(**data)
        else:
            raise TypeError(f"Invalid format for KeyValueEntry: '{data}'")

    def to_dict(self) -> dict:
        """Transform KeyValueEntry to dictionary.

        Creates a dictionary with the following keys: id, name, key, value, uri.

        Returns
        -------
        dict
        """
        return {
            "id": self.id,
            "name": self.name,
            "key": self.key,
            "value": self.value,
            "uri": self.uri,
        }

    def __str__(self) -> str:
        """Get string representation.

        Returns
        -------
        str
        """
        return str(self.to_dict())

    def __repr__(self) -> str:
        """Get string representation, including module and class name.

        Returns
        -------
        str
        """
        return (
            f"<{self.__class__.__module__}.{self.__class__.__qualname__}"
            f"({self.key}, {self.value}, {self.uri})>"
        )


class KeyValuePairs(UserDict):
    """A UserDict to store KeyValueEntries.

    Parameters
    ----------
    entries : Iterable
        an iterable containing entry information
    """

    def __init__(self, entries: Iterable[Union[Dict, KeyValueEntry]] = None):
        """Initialize the KeyValuePairs dictionary class.

        Parameters
        ----------
        entries: Iterable
            An iterable of dictionaries or KeyValueEntry, which will be inputted to the
            dictionary.
        """
        super().__init__()
        if entries:
            for item in entries:
                entry = KeyValueEntry.from_data(item)
                self.data[entry.key] = entry

    def __setitem__(self, key: str, item: Union[Dict, KeyValueEntry]) -> None:
        """Set item.

        Parameters
        ----------
        key: str
        item: dictionary or KeyValueEntry
        """
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
        """Get string representation, including module and class name.

        Returns
        -------
        str
        """
        return (
            f"{self.__class__.__module__}.{self.__class__.__qualname__}"
            f"({self.to_dict()!r})"
        )

    def _repr_html_(self) -> str:
        """Get HTML representation.

        Returns
        -------
        str
            HTML formatted string
        """
        return f"""<p><strong>KeyValuePairs</strong></p><p>{format_long_string(
            self.__str__(), 100)}</p>"""

    def to_dict(self) -> dict:
        """Get dictionary representation.

        Returns
        -------
        dict
            keys are the keys, and each value is the KeyValueEntry represented as
            a dict.
        """
        return {k: v.to_dict() for k, v in self.data.items()}
