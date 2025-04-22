"""The Metadata class that provides an interface to all types of cobra metadata."""

from collections import OrderedDict
from datetime import datetime
from typing import Dict, Iterable, List, Optional, Union

import cobra.core.metadata.custom as CA
import cobra.core.metadata.standardized as SA

from ..metadata.history import Creator, History


class Metadata:
    """Metadata of a cobrapy object.

    Metadata encodes additional information on an object, such as annotations. This
    information is mainly stored in the SBML annotation tags.

    Metadata consists of four components:
    - standardized: contains annotations in a standardized format, based on
      https://identifier.org identifiers and BioModels.net qualifiers.
      This information is also exposed via the object.annotation interface for
      full backwards compatibility with earlier cobrapy versions. See
      Metadata.standardized and Object.annotation.
    - history: storing the object history consisting of creators, created date, and
      modified dates.
    - custom: contains custom annotations as key-value pairs.
    - sbo - a single SBO term for the object.
    """

    def __init__(
        self,
        standardized: Optional[
            Union[
                Dict,
                List["SA.StandardizedAnnotation"],
                "SA.StandardizedAnnotationStore",
            ]
        ] = None,
        history: Optional[Union[Dict, History]] = None,
        custom: Optional[List] = None,
        sbo: str = "",
    ):
        """Initialize the Metadata class.

        Parameters
        ----------
        standardized : dict, list of StandardizedAnnotation, StandardizedAnnotationStore
            Collection of standardized annotations. Dictionaries and lists of
            StandardizedAnnotations are converted to a StandardizedAnnotationStore and
            set as the `Metadata.standardized` attribute.
        history : dict, History
            The history information as a History object or dictionary.
        custom: list
            Custom annotation key-value pairs.
        sbo: str
            The sbo term to use for the entity.
        """
        self.standardized = standardized
        self.history = history
        self.custom = CA.CustomAnnotationStore(custom)
        self.sbo = sbo

    @property
    def standardized(self) -> "SA.StandardizedAnnotationStore":
        """Get the standardized annotations.

        Returns
        -------
        StandardizedAnnotationStore
        """
        return self._standardized

    @standardized.setter
    def standardized(
        self,
        values: Optional[
            Union[
                Dict,
                Iterable["SA.StandardizedAnnotation"],
                "SA.StandardizedAnnotationStore",
            ]
        ],
    ) -> None:
        """Set the standardized annotations.

        Parameters
        ----------
        values: dict, list of StandardizedAnnoation or StandardizedAnnotationStore
            Lists and dicts are converted to StandardizedAnnotationStore using
            StandardizedAnnotationStore.from_data().

        See Also
        --------
        StandardizedAnnotationStore.from_data
        """
        self._standardized = SA.StandardizedAnnotationStore.from_data(values)

    def add_standardized(
        self, annotations: List[Union[Dict, "SA.StandardizedAnnotation"]]
    ) -> None:
        """Add one or more standardized annotations.

         This method will add StandardizedAnnotation objects to the standardized field.

         Parameters
         ----------
         annotations: list of dict or StandardizedAnnotation
            An iterable of CVTerm objects or CVTerm dict compatible objects.

        See Also
        --------
        StandardizedAnnotationStore.add
        """
        self._standardized.add(annotations)

    @property
    def history(self) -> History:
        """Get history of the object.

        Returns
        -------
        History
        """
        return self._history

    @history.setter
    def history(self, history: Optional[Union[Dict, History]]) -> None:
        """Set history of the object.

        Parameters
        ----------
        history: History or dict or None
            If None is given, will set the History to be empty.
            Dict is converted via History.from_data()
        """
        self._history = History.from_data(history)

    def add_creators(self, creators: Iterable[Union[Creator, dict]]) -> None:
        """Add one or more creators to the object history.

        The creators will be parsed and need to be a Creator object or a dictionary
        in the correct format.

        Parameters
        ----------
        creators:  Iterable of dict or Creator
            An iterable of dicts or Creator objects.

        See Also
        --------
        Creator.from_data()
        """
        self.history.creators.extend(
            [Creator.from_data(creator) for creator in creators]
        )

    def add_modification_dates(
        self, dates: Union[datetime, str, Iterable[Union[datetime, str]]]
    ) -> None:
        """Add modification dates to the object history.

        The dates will be parsed and need to be string in the acceptable format or
        datetime.

        Parameters
        ----------
        dates: Iterable or str or datetime
            An iterable of strings or datetime objects, or one str or one datetime.

        See Also
        --------
        History.parse_datetime()
        """
        if isinstance(dates, (str, datetime)):
            dates = [dates]
        self.history.modified_dates.extend([History.parse_datetime(d) for d in dates])

    @property
    def sbo(self) -> str:
        """Return the SBO term of the Metadata.

        Returns
        -------
        str: SBO as string
        """
        return self._sbo

    @sbo.setter
    def sbo(self, value: Union[str, List[str]]) -> None:
        """Set the SBO term."""
        if isinstance(value, list):
            value = value[0]
        self._sbo = value

    @property
    def custom(self) -> "CA.CustomAnnotationStore":
        """Returns the custom key-value pair annotations.

        Returns
        -------
        CustomAnnotationStore: The custom annotations.
        """
        return self._custom

    @custom.setter
    def custom(
        self,
        annotations: Union[
            Dict, List["CA.CustomAnnotation"], "CA.CustomAnnotationStore"
        ],
    ) -> None:
        """Set the custom key-value pair annotations.

        Parameters
        ----------
        annotations: dict, list of CustomAnnotation or CustomAnnotationStore
            A dictionary or CustomAnnotationStore instance that contains all custom
            annotation key-value pairs.
        """
        self._custom = CA.CustomAnnotationStore(annotations)

    # def __setitem__(self, key: str, value: Union[List, str]) -> None:
    #     """Set the item for accessing metadata as dict (the old style annotation).
    #
    #     Parameters
    #     ----------
    #     key: str
    #         provider key word.
    #     value: List or str
    #         A str that is one term or a list that will contain multiple terms
    #
    #     This function will first delete the existing value for the key, and then set
    #     it to the new value. Be careful - if you give this function incorrect input,
    #     the deletion will happen anyway, and the value of the key will be empty!
    #
    #     If the key is sbo, sets the self.sbo term to the first item in the list. The
    #     rest of the items in the list are ignored.
    #     Cobrapy support for multiple SBO terms is not implemented yet.
    #
    #     See Also
    #     --------
    #     `CVTermList().add_simple_annotations()`
    #     """
    #     if key == "sbo":
    #         if isinstance(value, list):
    #             value = value[0]
    #         self.sbo = value
    #     else:
    #         self._standardized[key] = value
    #
    # def __getitem__(
    #     self, key: Union[str, Qualifier]
    # ) -> Union[List[str], Dict[str, str]]:
    #     """Get item using old annotation type dictionary.
    #
    #     If the key is sbo, will return the sbo field directly.
    #     Otherwise, will query the annotations (old style) dictionary.
    #
    #     Note, that __setitem__, __getitem__ and __delitem__ will ignore custompairs.
    #     If you want to edit that field, use relevant functions for it.
    #
    #     Parameters
    #     ----------
    #     key: str
    #         provider key word.
    #
    #     Returns
    #     -------
    #     list
    #     """
    #     if key == "sbo" or key == "SBO":
    #         return [self.sbo]
    #     else:
    #         return self.standardized[key]
    #
    #
    # def __delitem__(self, key: str) -> None:
    #     """Delete item using old annotation type dictionary as reference.
    #
    #     If the key is sbo, will zero out the sbo field directly to an empty string.
    #     Otherwise, will delete from the annotations (old style) dictionary, using
    #     the value as a pattern to search in CVTermList resources.
    #     See `CVTermList.delete_simple_annotation()`
    #
    #     Parameters
    #     ----------
    #     key: str
    #         provider key word.
    #
    #     Note, that __setitem__, __getitem__ and __delitem__ will ignore custompairs.
    #     If you want to edit that field, use relevant functions for it.
    #     """
    #     if key == "sbo":
    #         self.sbo = ""
    #     else:
    #         self._standardized.delete_annotation(key)
    #
    def __eq__(self, other: Union[Dict, "Metadata"]) -> bool:
        """Compare two Metadata objects to find out whether they are equal.

        If given a dict, the dictionary is converted to Metadata and then compared.

        Two metadata objects are equal (the function will return True) if
        - standardized annotatios are equal
        - all attributes of the history are equal
        - custom annotations are equal
        If one of these three conditions is not true, the function will return False.

        Parameters
        ----------
        other: Metadata or dict

        Returns
        -------
        bool: True if equal, False otherwise.
        """
        if isinstance(other, dict):
            return self == Metadata.from_dict(other)
        elif isinstance(other, Metadata):
            return (
                (self.standardized == other.standardized)
                and (self.history == other.history)
                and (self.custom == other.custom)
            )
        else:
            raise TypeError(
                "Can only compare Metadata objects to dictionaries or other Metadata"
                f"objects, not: {type(other)}."
            )

    def __ne__(self, other) -> bool:
        """Compare two Metadata objects to find out whether they are not equal.

        Returns the inverse of `Metatdata.__eq__`.

        Parameters
        ----------
        other: Metadata or dict

        Returns
        -------
        bool: False if equal, True otherwise.

        See Also
        --------
        Metadata.__eq__()
        """
        return not self.__eq__(other)

    #
    # def __iter__(self) -> Iterator:
    #     """Iterate over the Metadata annotations dict.
    #
    #     This function will iterate over the annotations (old-style) dictionary.
    #
    #     Returns
    #     -------
    #     Iterator
    #     """
    #     return iter(self.annotations)
    #
    # def __len__(self) -> int:
    #     """Return the length of the Metadata annotations dict.
    #
    #     The length of the annotations dict will be returned as a int.
    #
    #     Returns
    #     -------
    #     int
    #         result of running len on the annotations dict
    #     """
    #     return len(self.annotations)
    #
    # def __str__(self) -> str:
    #     """Return the Metadata as str.
    #
    #     The annotations dict will be returned as a string. This does not include all
    #     of the possible Metadata fields.
    #
    #     Returns
    #     -------
    #     str
    #     """
    #     return str(dict(self.annotations))
    #
    # def __repr__(self) -> str:
    #     """Return the Metadata as str with module, class, and code to recreate it.
    #
    #     If Metadata has standardized, history or custompairs, the dictionary will
    #     contain them as keys.
    #     If not, the dictionary will have only the annotations dictionary.
    #     In either case, if sbo field is set, it will be outputted in the dictionary.
    #
    #     Returns
    #     -------
    #     str
    #
    #     See Also
    #     --------
    #     Metadata.from_dict()
    #     """
    #     repr_str = (
    #         f"{self.__class__.__module__}.{self.__class__.__qualname__}.from_dict("
    #     )
    #     cvterms = self.standardized
    #     history = self.history
    #     custompairs = self.custompairs
    #
    #     if cvterms or not history.is_empty() or custompairs:
    #         repr_str += (
    #             f"'standardized': {self.standardized.to_list_of_dicts()},"
    #             f"'history': {self.history.to_dict()},"
    #             f"'custompairs': {self.custompairs},"
    #         )
    #         if self.sbo:
    #             repr_str += f"'sbo': {self.sbo})"
    #     else:
    #         anno_dict = self.annotations
    #         if self.sbo:
    #             anno_dict.update({"sbo": self.sbo})
    #         repr_str += f"{anno_dict})"
    #     return repr_str
    #
    def to_dict(self) -> Dict:
        """Create a dictionary from the Metadata object.

        The dictionary will contain any of the keys 'sbo', 'standardized', 'history',
        and 'custom', if the corresponding attributes are not empty.

        Returns
        -------
        dict

        See Also
        --------
        StandardizedAnnotationStore.to_list_of_dicts
        History.to_dict
        CustomAnnotationStore.to_dict
        """
        d = OrderedDict()
        if self.sbo:
            # set first SBO term as sbo
            d["sbo"] = self.sbo

        if self.standardized:
            d["standardized"] = self.standardized.to_list_of_dicts()

        if self.history and not self.history.is_empty():
            d["history"] = self.history.to_dict()

        if self.custom:
            d["custom"] = self.custom.to_dict()

        return d

    @staticmethod
    def from_dict(data: Dict) -> "Metadata":
        """Generate a Metadata instance from dictionary.

        The dictionary should have any of the keys 'standardized', 'history', 'sbo', and
        'custompairs', which will be converted to the corresponding attributes.

        Parameters
        ----------
        data: dict
            Dictionary to transform into Metadata.

        Returns
        -------
        Metadata
        """
        standardized = data.get("standardized", None)
        history = data.get("history", None)
        custom = data.get("custom", None)

        if standardized or history or custom:
            annotation = Metadata(
                standardized=standardized, history=history, custom=custom
            )
        else:
            annotation = Metadata()
            # annotation.standardized.add_simple_annotations(data)
            # raise ValueError()
            # TODO: Fix

        if "sbo" in data:
            annotation.sbo = data["sbo"]

        return annotation
