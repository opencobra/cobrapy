from collections import OrderedDict, MutableMapping
from typing import Dict, Iterator, List, Union

from cobra.core.metadata.cvterm import CVTerms
from cobra.core.metadata.history import History, Creator, HistoryDatetime
from cobra.core.metadata.keyvaluepairs import KeyValuePairs


class MetaData(MutableMapping):
    """Meta-data of a cobrapy object.

    Meta-data encodes additional information on an object such as annotations
    or notes. Such information is currently stored in SBML on the notes and
    annotation fields.

    Meta-data consists of three components:
    - CVTerms: storing resource:identifier annotation information. The annotation
      information is exposed via the dict interface for full backwards compatibility
      to the earlier object.annotation field.
    - History: storing the object history consisting of creators, created date, and
               modified dates.
    - KeyValuePairs, a list of key-value pairs to store additional information

    Parameters
    ----------
    cvterms : dict or CVTerms object
        The cvterms store annotations to external resources
    history : dict, History
        The history stores information about the creator,
        created and modified dates.
    keyvaluepairs : list
        Key-value pairs which are not suitable to be
        represented anywhere else in the model.
        Data is represented as an OrderedDict.
    """

    def __init__(
        self,
        cvterms: Union[Dict, CVTerms] = None,
        history: Union[Dict, History] = None,
        keyvaluepairs: List = None,
    ):
        self._cvterms = None  # type: CVTerms
        self._history = None  # type: History
        self._keyvaluepairs = None  # type: KeyValuePairs

        # use setters
        self.cvterms = cvterms
        self.history = history
        self.keyvaluepairs = keyvaluepairs

    @property
    def annotations(self) -> Dict:
        """Backwards compatible annotations."""
        return self.cvterms.annotations

    def __setitem__(self, key: str, value: List) -> None:
        self._cvterms.add_simple_annotations(dict({key: value}))

    def __getitem__(self, key: str) -> None:
        return self.annotations[key]

    def __delitem__(self, key: str) -> None:
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
    def cvterms(self) -> "CVTerms":
        return self._cvterms

    @cvterms.setter
    def cvterms(self, cvterms: Union[Dict, CVTerms]) -> None:
        self._cvterms = CVTerms.from_data(cvterms)

    def add_cvterms(self, cvterms: Union[Dict, CVTerms]) -> None:
        self._cvterms.add_cvterms(cvterms)

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
        return self._keyvaluepairs

    @keyvaluepairs.setter
    def keyvaluepairs(self, keyvaluepairs: Union[Dict, KeyValuePairs]) -> None:
        self._keyvaluepairs = KeyValuePairs(keyvaluepairs)

    def to_dict(self) -> Dict:
        """Creates string dictionary for serialization"""
        d = OrderedDict()
        if "sbo" in self and self["sbo"] != []:
            # set first SBO term as sbo
            d["sbo"] = self["sbo"][0]

        if self.cvterms:
            d["cvterms"] = self.cvterms.to_dict()

        if self.history and not self.history.is_empty():
            d["history"] = self.history.to_dict()

        if self.keyvaluepairs:
            d["keyvaluepairs"] = self.keyvaluepairs.to_dict()

        return d
