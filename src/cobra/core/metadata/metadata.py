"""Define the cobra MetaData class."""

from collections import OrderedDict
from collections.abc import MutableMapping
from datetime import datetime
from typing import Dict, Iterable, Iterator, List, Optional, Union

from ..metadata.custompairs import KeyValuePairs
from ..metadata.cvterm import CVTerm, CVTermList
from ..metadata.history import Creator, History


class MetaData(MutableMapping):
    """Meta-data of a cobrapy object.

    Meta-data encodes additional information on an object such as annotations
    or notes. Such information is currently stored in SBML on the notes and
    annotation fields.

    Meta-data consists of four components:
    - standardized: storing resource:identifier annotation information in a CVTermList
      object. The annotation information is exposed via the dict interface for
      full backwards compatibility to the earlier object.annotation field. See
      CVTermList.annotations() and MetaData.annotations()
    - sbo - a single SBO term for the object. Can theoretically be placed in
      standardized, but the most descriptive term should be kept in sbo. Right now,
      cobrapy doesn't allow multiple sbo terms.
    - History: storing the object history consisting of creators, created date, and
               modified dates.
    - custompairs (KeyValuePairs), a list of key-value pairs to
      store additional information

    Parameters
    ----------
    cvterms : dict or CVTermList object
        The standardized store annotations to external resources
    history : dict, History
        The history stores information about the creator,
        created and modified dates.
    sbo: str
        The sbo term to use for the entity. If you want to use more than one SBO term
        (not recommended), use SBO in identifers.org format and put it in standardized.
    keyvaluepairs : list
        Key-value pairs which are not suitable to be
        represented anywhere else in the model.
        Data is represented as an OrderedDict.
    """

    def __init__(
        self,
        cvterms: Union[Dict, CVTermList] = None,
        history: Union[Dict, History] = None,
        sbo: str = "",
        keyvaluepairs: List = None,
    ):
        """Initialize the MetaData class.

        Parameters
        ---------
        cvterms : dict or CVTermList, optional
            Which controlled vocabulary terms does the metadata have. Default None.
        history: dict or History, optional
            History of annotation, including creation data, creators, and optional
            modificiation date. Default None.
        sbo: str
            SBO term, if relevant. Default "".
        keyvaluepairs: KeyValuePairs
            For annotations that don't match the identifiers.org format.

        """
        self._standardized = CVTermList.from_data(cvterms)
        self._history = History.from_data(history)
        self._custompairs = KeyValuePairs(keyvaluepairs)
        self._sbo = sbo

        # use setters
        self.sbo = sbo
        self.standardized = cvterms
        self.history = history
        self.custompairs = keyvaluepairs

    @property
    def standardized(self) -> "CVTermList":
        """Get standardized field of MetaData.

        Returns
        -------
        CVTermList
        """
        return self._standardized

    @standardized.setter
    def standardized(self, cvterms: Union[Dict, CVTermList]) -> None:
        """Set standardized field of MetaData with controlled vocabulary (CVTerm).

        Parameters
        ----------
        cvterms: dict or CVTermList
            dict is converted to CVTermList using CVTermList.from_data().
            The wrong type will lead to a TypeError being raised.
        """
        self._standardized = CVTermList.from_data(cvterms)

    def add_cvterms(self, cvterms: Iterable[Union[Dict, CVTerm]]) -> None:
        """Add one or more CVTerm objects to the standardized field.

         This method will add CVTerm objects to the standardized field.

         Parameters
         ----------
         cvterms: Iterable
            An iterable of CVTerm objects or CVTerm dict compatible objects.
            CVTermList is an acceptable Iterable in this case.

        See Also
        --------
        CVTermList.add_cvterms()
        """
        self._standardized.add_cvterms(cvterms)

    @property
    def annotations(self) -> Dict:
        """Backwards compatible annotations."""
        anno_dict = self.standardized.annotations
        if self.sbo:
            anno_dict["sbo"] = self.sbo
        return anno_dict

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
        """Return the SBO term of the MetaData.

        Returns
        -------
        str: SBO as string
        """
        return self._sbo

    @sbo.setter
    def sbo(self, value: str) -> None:
        """Set the SBO term."""
        if isinstance(value, list):
            value = value[0]
        self._sbo = value

    @property
    def keyvaluepairs(self) -> KeyValuePairs:
        return self._custompairs

    @keyvaluepairs.setter
    def keyvaluepairs(self, keyvaluepairs: Union[Dict, KeyValuePairs]) -> None:
        self._custompairs = KeyValuePairs(keyvaluepairs)

    def __setitem__(self, key: str, value: Union[List, str]) -> None:
        """Set the item for accessing metadata as dict (the old style annotation).

        Parameters
        ----------
        key: str
            provider key word.
        value: List or str
            A str that is one term or a list that will contain multiple terms

        This function will first delete the existing value for the key, and then set
        it to the new value. Be careful - if you give this function incorrect input,
        the deletion will happen anyway, and the value of the key will be empty!

        If the key is sbo, sets the self.sbo term to the first item in the list. The
        rest of the items in the list are ignored.
        Cobrapy support for multiple SBO terms is not implemented yet.

        See Also
        --------
        `CVTermList().add_simple_annotations()`
        """
        if key == "sbo":
            if isinstance(value, list):
                value = value[0]
            self.sbo = value
        else:
            self.__delitem__(key)
            self._standardized.add_simple_annotations(dict({key: value}))

    def __getitem__(self, key: str) -> List:
        """Get item using old annotation type dictionary.

        If the key is sbo, will return the sbo field directly.
        Otherwise, will query the annotations (old style) dictionary.

        Note, that __setitem__, __getitem__ and __delitem__ will ignore custompairs. If
        you want to edit that field, use relevant functions for it.

        Parameters
        ----------
        key: str
            provider key word.

        Returns
        -------
        list
        """
        if key == "sbo":
            return [self.sbo]
        else:
            return self.annotations[key]

    def __delitem__(self, key: str) -> None:
        """Delete item using old annotation type dictionary as reference.

        If the key is sbo, will zero out the sbo field directly to an empty string.
        Otherwise, will delete from the annotations (old style) dictionary, using
        the value as a pattern to search in CVTermList resources.
        See `CVTermList.delete_simple_annotation()`

        Parameters
        ----------
        key: str
            provider key word.

        Note, that __setitem__, __getitem__ and __delitem__ will ignore custompairs. If
        you want to edit that field, use relevant functions for it.
        """
        if key == "sbo":
            self.sbo = ""
        else:
            self._standardized.delete_annotation(key)

    def __eq__(self, other: Union["MetaData", Dict]) -> bool:
        """Compare two MetaData objects to find out whether they are the same.

        If given a dict, the dictionary is converted to MetaData and then compared.

        Equality is defined in two ways, depending if history and custompairs exist.
        1) If history or custompairs exist and are not empty/None, the two objects are
           equal (the function will return True) if
           - standardized (CVTermList) are equal
           - all attributes of the history (History) are equal
           - custompairs are identical
           If one of these three conditions is not true, the function will return False.
        2) If only standardized exist, while history and custompairs are empty
           (history.is_empty() is True and custompairs is None) for both objects,
           then the annotations field (CVTermList.annotations()) is compared as two
           dictionaries.

        Parameters
        ----------
        other: MetaData or dict

        Returns
        -------
        bool: True if equal, False otherwise.
        """
        if isinstance(other, dict):
            return self == MetaData.from_dict(other)
        if (
            self.standardized
            and self.history.is_empty()
            and not self.custompairs
            and other.standardized
            and other.history.is_empty()
            and not other.custompairs
        ):
            return self.annotations == other.annotations

        return (
            (self.standardized == other.standardized)
            and (self.history == other.history)
            and (self.custompairs == other.custompairs)
        )

    def __ne__(self, other) -> bool:
        """Define not-equal to override default.

        Parameters
        ----------
        other: MetaData or dict

        Returns
        -------
        bool: False if equal, True otherwise.

        See Also
        --------
        MetaData.__eq__()
        """
        return not self.__eq__(other)

    def __iter__(self) -> Iterator:
        """Iterate over the MetaData annotations dict.

        This function will iterate over the annotations (old-style) dictionary.

        Returns
        -------
        Iterator
        """
        return iter(self.annotations)

    def __len__(self) -> int:
        """Return the length of the MetaData annotations dict.

        The length of the annotations dict will be returned as a int.

        Returns
        -------
        int
            result of running len on the annotations dict
        """
        return len(self.annotations)

    def __str__(self) -> str:
        """Return the MetaData as str.

        The annotations dict will be returned as a string. This does not include all
        of the possible MetaData fields.

        Returns
        -------
        str
        """
        return str(dict(self.annotations))

    def __repr__(self) -> str:
        """Return the MetaData as str with module, class, and code to recreate it.

        If MetaData has standardized, history or custompairs, the dictionary will
        contain them as keys.
        If not, the dictionary will have only the annotations dictionary.
        In either case, if sbo field is set, it will be outputted in the dictionary.

        Returns
        -------
        str

        See Also
        --------
        MetaData.from_dict()
        """
        repr_str = (
            f"{self.__class__.__module__}.{self.__class__.__qualname__}.from_dict("
        )
        cvterms = self.standardized
        history = self.history
        keyValuepairs = self.keyvaluepairs

        if cvterms or not history.is_empty() or keyValuepairs:
            repr_str += (
                f"'standardized': {self.standardized.to_list_of_dicts()},"
                f"'history': {self.history.to_dict()},"
                f"'custompairs': {self.custompairs},"
            )
            if self.sbo:
                repr_str += f"'sbo': {self.sbo})"
        else:
            anno_dict = self.annotations
            if self.sbo:
                anno_dict.update({"sbo": self.sbo})
            repr_str += f"{anno_dict})"
        return repr_str

    def to_dict(self) -> Dict:
        """Create string dictionary for serialization.

        Returns
        -------
        dict
        """
        d = OrderedDict()
        if self.sbo:
            # set first SBO term as sbo
            d["sbo"] = self.sbo

        if self.standardized:
            d["standardized"] = self.standardized.to_list_of_dicts()

        if self.history and not self.history.is_empty():
            d["history"] = self.history.to_dict()

        if self.keyvaluepairs:
            d["custompairs"] = self.keyvaluepairs.to_dict()

        return d

    @staticmethod
    def from_dict(data: Dict) -> "MetaData":
        """Generate MetaData from dictionary.

        This function has two options
        1) If the dictionary has standardized, history, or custompairs, it uses
           the relevant fields to create MetaData.
        2) If the dictionary has none of the above fields, it creates an empty
           MetaData object, and then fills it using CVTermList.add_simple_annotation()
           where the keys are namespaces and the values are identifiers.
           In this case, all qualifiers will be "bqb_is".

        In either case, if "sbo" is present as one of the keys, this will fill the sbo
        field of the MetaData using the value of "sbo" key.

        Parameters
        ----------
        data: dict
            Dictionary to transform into MetaData.

        Returns
        -------
        MetaData
        """
        cvterms = data.get("standardized", None)
        history = data.get("history", None)
        keyValuepairs = data.get("custompairs", None)

        if cvterms or history or keyValuepairs:
            annotation = MetaData(cvterms, history, keyValuepairs)
        else:
            annotation = MetaData()
            annotation.standardized.add_simple_annotations(data)

        if "sbo" in data:
            annotation["sbo"] = data["sbo"]

        return annotation
