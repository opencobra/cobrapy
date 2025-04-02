"""Class to deal with Key-Value pairs.

Key-Value pairs are described in the SBML FBC3 proposal. For the latest
version of the FBC3 proposal, see Release Candidate 1:
https://github.com/sbmlteam/sbml-specifications/blob/develop/sbml-level-3/version-1/fbc/spec/sbml-fbc-version-3-release-1.pdf
"""

# TODO: Update docstring with final release, when available.
import uuid
from collections import UserDict
from dataclasses import asdict
from typing import Dict, Iterable, Optional, Union

from ...util import format_long_string
from .. import object as cobject


class CustomAnnotation(cobject.Object):
    """Single key-value entry.

    The key is an attribute on the entry.

    Parameters
    ----------
    key: str
        Defined as mandatory in the FBC3 standard.
    value: str
        optional. Default None.
    uri: str
        Can be a URN or URL. Optional (default None).
    """

    def __init__(
        self, key: str, value: Optional[str] = None, uri: Optional[str] = None
    ):
        super(__class__, self).__init__()
        self._key = key
        self._value = value
        self._uri = uri
        self._target = None

    def _set_target(self, target: Optional["cobject.Object"]) -> None:
        self._target = target

    def remove_from_object(self):
        if self._target is None:
            raise ValueError(
                "Cannot remove annontation, since no object is associated with annotation."
            )
        self._target.remove_annotations(self)

    # We should probably make key read-only, to prevent keys from becoming duplicate in
    # CustomAnnotationList objects.
    @property
    def key(self) -> str:
        return self._key

    @property
    def value(self) -> Optional[str]:
        return self._value

    @value.setter
    def value(self, value: Optional[str]) -> None:
        self._value = value

    @property
    def uri(self) -> Optional[str]:
        return self._uri

    @uri.setter
    def uri(self, uri: Optional[str]) -> None:
        self._uri = uri

    @staticmethod
    def from_data(
        data: Optional[Union[Dict, "CustomAnnotation"]],
    ) -> "CustomAnnotation":
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
        if isinstance(data, CustomAnnotation):
            return data
        elif isinstance(data, dict):
            if "key" not in data:
                data["key"] = uuid.uuid4().hex
            ann = CustomAnnotation(
                key=data["key"], value=data.get("value"), uri=data.get("uri")
            )
            if "id" in data:
                ann.id = data["id"]
            if "name" in data:
                ann.name = data["name"]
            # TODO: Handle annotations
            return ann
        else:
            raise TypeError(f"Invalid format for CustomAnnotation: '{data}'")

    def _asdict(self) -> dict:
        return {
            k: v
            for k in ["key", "value", "uri", "id", "name"]
            if (v := getattr(self, k, None)) is not None and v != ""
        }

    def __str__(self) -> str:
        """Get string representation of the KeyValueEntry as dictionary.

        Returns
        -------
        str
        """
        return str(self._asdict())

    def __repr__(self) -> str:
        """Get string representation, including module and class name.

        Returns
        -------
        str
        """
        return (
            f"{self.__class__.__module__}.{self.__class__.__qualname__}"
            f"({repr(self.key)}, {repr(self.id)}, {repr(self.name)}, {repr(self.value)}"
            f", {repr(self.uri)})"
        )


class CustomAnnotationList(UserDict):
    """A UserDict to store KeyValueEntries.

    Parameters
    ----------
    entries : Iterable
        an iterable containing entry information
    """

    def __init__(
        self,
        entries: Optional[
            Union[Iterable[Union[Dict, CustomAnnotation]], "CustomAnnotationList"]
        ] = None,
    ):
        """Initialize the KeyValuePairs dictionary class.

        Parameters
        ----------
        entries: Iterable
            An iterable of dictionaries or KeyValueEntry, which will be inputted to the
            dictionary.
        """
        super().__init__()
        if entries is None:
            return
        elif isinstance(entries, CustomAnnotationList):
            self.data = entries.data.copy()
        else:
            for item in entries:
                entry = CustomAnnotation.from_data(item)
                self.data[entry.key] = entry

    def __setitem__(self, key: str, item: Union[Dict, CustomAnnotation]) -> None:
        """Set item.

        Parameters
        ----------
        key: str
        item: dictionary or KeyValueEntry
        """
        entry = CustomAnnotation.from_data(item)
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
        return {k: asdict(v) for k, v in self.data.items()}

    def add(
        self,
        items: Union[CustomAnnotation, Dict, Iterable[Union[CustomAnnotation, Dict]]],
    ) -> None:
        if isinstance(items, CustomAnnotation):
            items = [items]
        elif isinstance(items, dict):
            items = [CustomAnnotation.from_data(items)]
        for item in items:
            item = CustomAnnotation.from_data(item)
            self.data[item.key] = item

    def remove(
        self,
        items: Union[CustomAnnotation, str, Iterable[Union[CustomAnnotation, str]]],
    ) -> None:
        if isinstance(items, CustomAnnotation):
            items = [items.key]
        elif isinstance(items, str):
            items = [items]

        for item in items:
            if isinstance(item, CustomAnnotation):
                item = item.key
            # If CustomAnnotation object is removed from CustomAnnotationList, it will
            # also not belong to the target Object anymore.
            self.data[item]._set_target(None)
            del self.data[item]

    # query

    # add_key_value_pair
    # delete_key_value_pair?? Maybe with query
