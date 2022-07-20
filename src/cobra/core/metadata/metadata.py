"""Define the cobra MetaData class."""

from collections import OrderedDict
from collections.abc import MutableMapping
from typing import Dict, Iterable, Iterator, List, Union

from ..metadata.cvterm import CVTerm, CVTermList
from ..metadata.history import Creator, History
from ..metadata.custompairs import KeyValuePairs


class MetaData(MutableMapping):
    """Meta-data of a cobrapy object.

    Meta-data encodes additional information on an object such as annotations
    or notes. Such information is currently stored in SBML on the notes and
    annotation fields.

    Meta-data consists of three components:
    - CVTermList: storing resource:identifier annotation information. The annotation
      information is exposed via the dict interface for full backwards compatibility
      to the earlier object.annotation field.
    - History: storing the object history consisting of creators, created date, and
               modified dates.
    - KeyValuePairs, a list of key-value pairs to store additional information

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
            If there are ????
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
    def annotations(self) -> Dict:
        """Backwards compatible annotations."""
        anno_dict = self.standardized.annotations
        if self.sbo:
            anno_dict["sbo"] = self.sbo
        return anno_dict

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
        self._sbo = value

    def __setitem__(self, key: str, value: Union[List, str]) -> None:
        """Set the item for accessing metadata as dict (the old style annotation).

        Parameters
        ----------
        key: str
            provider key word.
        value: List
            A list that will contain either term(s) or qualifier and term(s).

        If the key is sbo, sets the self.sbo term to the first item in the list. If
        you'd like to add multiple SBO terms, use the CVTermList() and add sbo as
        identifiers.org formatted external resources/links.

        See Also
        --------
        CVTermList().add_simple_annotations()

        """
        if key == "sbo":
            if isinstance(value, list):
                value = value[0]
            self.sbo = value
        else:
            self._standardized.add_simple_annotations(dict({key: value}))

    def __getitem__(self, key: str) -> Union[str, List]:
        if key == "sbo":
            return [self.sbo]
        else:
            return self.annotations[key]

    def __delitem__(self, key: str) -> None:
        if key == "sbo":
            self.sbo = ""
        else:
            del self.annotations[key]

    def __iter__(self) -> Iterator:
        return iter(self.annotations)

    def __len__(self) -> int:
        return len(self.annotations)

    def __str__(self) -> str:
        return str(dict(self.annotations))

    def __repr__(self) -> str:
        return str(dict(self.annotations))

    @property
    def standardized(self) -> "CVTermList":
        return self._standardized

    @standardized.setter
    def standardized(self, cvterms: Union[Dict, CVTermList]) -> None:
        self._standardized = CVTermList.from_data(cvterms)

    def add_cvterms(self, cvterms: Iterable[Union[Dict, CVTerm]]) -> None:
        self._standardized.add_cvterms(cvterms)

    @property
    def history(self) -> History:
        return self._history

    @history.setter
    def history(self, history: Union[Dict, History]) -> None:
        self._history = History.from_data(history)

    def add_creator(self, creator: Creator):
        self.history.creators.append(creator)

    @property
    def keyvaluepairs(self) -> KeyValuePairs:
        return self._custompairs

    @keyvaluepairs.setter
    def keyvaluepairs(self, keyvaluepairs: Union[Dict, KeyValuePairs]) -> None:
        self._custompairs = KeyValuePairs(keyvaluepairs)

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
        cvterms = data.get("standardized", None)
        history = data.get("history", None)
        keyValuepairs = data.get("custompairs",  None)

        if cvterms or history or keyValuepairs:
            annotation = MetaData(cvterms, history, keyValuepairs)
        else:
            annotation = MetaData()
            annotation.standardized.add_simple_annotations(data)

        if "sbo" in data:
            annotation["sbo"] = data["sbo"]

        return annotation
