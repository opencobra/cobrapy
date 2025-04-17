"""Classes to handle custom annotations.

Custom annotations correspond to key-value pairs, described in the SBML FBC3 proposal.
For the latest version of the FBC3 proposal, see Release Candidate 1:
https://github.com/sbmlteam/sbml-specifications/blob/develop/sbml-level-3/version-1/fbc/spec/sbml-fbc-version-3-release-1.pdf
"""

# TODO: Update docstring with final release, when available.
import uuid
from collections import UserDict
from typing import Dict, Iterable, Optional, Union

from ...util import format_long_string
from .. import object as cobject


class CustomAnnotation(cobject.Object):
    """Custom annotation entry that represents a key-value/URI pair.

    `CustomAnnotation` object allow users to associate an arbitrary string value and/or
    URI with a key. This key-value store is meant to supplement the
    StandardizedAnnotation class, where MIRIAM-compliant annotations can be associated
    with model components. Whenever possible, `StandardizedAnnotation` objects should be
    preferred over CustomAnnotation objects, since they enable better interoperability
    between modelling tools.
    """

    def __init__(
        self, key: str, value: Optional[str] = None, uri: Optional[str] = None
    ):
        """Initialize a `CustomAnnotation` key-value pair object.

         Parameters
        ----------
        key: str
            Mandatory and descriptive string used for accessing the value/URI.
        value: str, optional
            Associated value. Default None.
        uri: str, optional
            URN or URL that links to the definition of the used key. Default None.
        """
        super(__class__, self).__init__()
        self._key = key
        self._value = value
        self._uri = uri
        self._parent = None

    def _set_parent(self, parent: Optional["CustomAnnotationStore"]) -> None:
        self._parent = parent

    def remove_from_object(self):
        """Remove this CustomAnnotation object from its `CustomAnnotationStore`.

        This method only removes and disassociates this object from its associated
        `CustomAnnotationStore` (e.g. `component.metadata.custom`), it does not
        delete the object itself.

        Raises
        ------
        ValueError if the object is not associated with a `CustomAnnotationStore`.

        See Also
        --------
        CustomAnnotationStore.add
        Object.add_annotation
        """
        if self._parent is None:
            raise ValueError(
                "Cannot remove annotation since no object is associated with it."
            )
        self._parent.remove(self)

    @property
    def key(self) -> str:
        """Get the key.

        The key is read-only to prevent duplication of keys in `CustomAnnotationStore`
        objects.

        Returns
        -------
        str
            The key of the `CustomAnnotation` object.
        """
        return self._key

    @property
    def value(self) -> Optional[str]:
        """Get the value associated to the key.

        Returns
        -------
        str or None
            The value associated with the key, or None if no value was provided.
        """
        return self._value

    @value.setter
    def value(self, value: Optional[str]) -> None:
        """Set the value associated to the key.

        Parameters
        ----------
        value: str or None
            The value to be associated with the key.
        """
        self._value = value

    @property
    def uri(self) -> Optional[str]:
        """Get the URI (URN or URL) associated to the key.

        Returns
        -------
        str or None
            The URI associated with the key, or None if no URI was provided.
        """
        return self._uri

    @uri.setter
    def uri(self, uri: Optional[str]) -> None:
        """Set the URI associated to the key.

        Parameters
        ----------
        uri: str or None
            The URI to be associated with the key.
        """
        self._uri = uri

    @staticmethod
    def from_data(
        data: Union[Dict, "CustomAnnotation"],
    ) -> "CustomAnnotation":
        """Create a `CustomAnnotation` object from data.

        Parameters
        ----------
        data: dict or CustomAnnotation
            Data to use to create the `CustomAnnotation`object. If data is of type dict,
            it should contain the key "key" and optionally "value", "uri", "id" and
            "name". If the data is already a `CustomAnnotation` object, this object will
            be returned.

        Returns
        -------
        CustomAnnotation

        Raises
        ------
        TypeError if data is not a dict or CustomAnnotation object.

        See Also
        --------
        to_dict
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

    def to_dict(self) -> dict:
        """Create a dictionary with the data of the `CustomAnnotation` object.

        Returns
        -------
        dict
            Dictionary containing all the `CustomAnnotation` data. Dictionary will
            contain the key "key" and optionally "value", "uri", "id" and "name".
        """
        return {
            k: v
            for k in ["key", "value", "uri", "id", "name"]
            if (v := getattr(self, k, None)) is not None and v != ""
        }

    def __str__(self) -> str:
        """Get string representation of the CustomAnnotation as dictionary.

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
            f"{self.__class__.__module__}.{self.__class__.__qualname__}"
            f"({repr(self.key)}, {repr(self.id)}, {repr(self.name)}, {repr(self.value)}"
            f", {repr(self.uri)})"
        )


class CustomAnnotationStore(UserDict):
    """A dict-like object that stores a collection of `CustomAnnotation` objects.

    This store allows users to associate arbitrary string values and/or
    URIs with keys. This key-value store is meant to supplement the standardized
    MIRIAM-compliant annotations. Whenever possible, standardized annotations are
    preferred over custom annotations, since they enable better interoperability
    between modelling tools.
    """

    def __init__(
        self,
        entries: Optional[
            Union[Iterable[Union[Dict, CustomAnnotation]], "CustomAnnotationStore"]
        ] = None,
    ):
        """Initialize the dict-like CustomAnnotationStore class.

        Parameters
        ----------
        entries: None, CustomAnnontationStore or list of CustomAnnotation or dicts,
        optional
            Custom annotations to initialize the store with. Default None.
        """
        super().__init__()
        if entries is None:
            return
        elif isinstance(entries, CustomAnnotationStore):
            self.data = entries.data.copy()
        else:
            for item in entries:
                entry = CustomAnnotation.from_data(item)
                self.data[entry.key] = entry

    def __setitem__(
        self, key: str, item: Optional[Union[Dict, CustomAnnotation, str]]
    ) -> None:
        """Set the value and/or URI associated with the key.

        Parameters
        ----------
        key: str
            Key used to look up value/URI in store.
        item: dict, str, CustomAnnotation or None
            Value/URI to associate to key. If item is of type str, this will be
            interpreted as value. If a dict is provided, this dict can contain any of
            the keys "value", "uri", "id", "name", "key".

        Raises
        ------
        ValueError if the "key" key in the provided dictionary or the key attribute of
            the CustomAnnotation does not match the key argument.

        See Also
        --------
        add
        """
        if isinstance(item, dict):
            if "key" in item and item["key"] != key:
                raise ValueError(
                    "The key in the annotation dictionary is not equal to the key "
                    "provided in the index."
                )
            # Make sure the key is also provided to the CustomAnnotation class.
            item = item | {"key": key}
        elif isinstance(item, CustomAnnotation):
            if item.key != key:
                raise ValueError(
                    "The key in the annotation object is not equal to the key "
                    "provided in the index."
                )
        elif isinstance(item, str):
            item = {"key": key, "value": item}
        else:
            raise TypeError(
                "CustomAnnotationStore entries should be provided as dict, str "
                "or CustomAnnotation object."
            )
        self.add(item, overwrite=True)

    def __str__(self) -> str:
        """Convert the custom annotation store to str.

        Returns
        ------
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

    def add(
        self,
        items: Union[CustomAnnotation, Dict, Iterable[Union[CustomAnnotation, Dict]]],
        overwrite: bool = False,
    ) -> None:
        """Add custom annotations to the store.

        Parameters
        ----------
        items: dict, CustomAnnotation or a list of dict or CustomAnnotation objects
            Custom annotation items to add to the store. CustomAnnotation objects can
            either be provided directly, or as a dictionary.
        overwrite: bool, optional
            Whether to overwrite an existing custom annotation with the same key.
            Default False.

        Raises
        ------
        IndexError if overwrite is False and the key of a custom annotation already
            exists in the store.
        """
        if isinstance(items, CustomAnnotation):
            items = [items]
        elif isinstance(items, dict):
            items = [CustomAnnotation.from_data(items)]
        for item in items:
            item = CustomAnnotation.from_data(item)
            if not overwrite and item.key in self.data:
                raise IndexError(f"Key '{item.key}' already exists in store.")
            self.data[item.key] = item

    def remove(
        self,
        items: Union[CustomAnnotation, str, Iterable[Union[CustomAnnotation, str]]],
    ) -> None:
        """Remove a custom annotation from the store.

        Parameters
        ----------
        items: str, CustomAnnotation or a list of str or CustomAnnotation objects
            Remove the proved CustomAnnotation objects from the store. If str objects
            are provided, they are interpreted as keys of the annotations to remove.

        Raises
        ------
        ValueError if a provided CustomAnnotation object is not present in the store.
        """
        if isinstance(items, (str, CustomAnnotation)):
            items = [items]

        for item in items:
            if isinstance(item, CustomAnnotation):
                ann = self.data[item.key]
                if ann is not item:
                    raise ValueError(
                        "Provided custom annotation does not match the "
                        "custom annotation in the store."
                    )
                item = ann
            # If CustomAnnotation object is removed from CustomAnnotationStore, it will
            # also not belong to the parent object anymore.
            self.data[item]._set_parent(None)
            del self.data[item]

    # query

    # add_key_value_pair
    # delete_key_value_pair?? Maybe with query
