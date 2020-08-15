import collections
from collections import defaultdict
from typing import Dict, Iterator, List, Union

from cobra.core.metadata.cvterm import CVTerms
from cobra.core.metadata.history import History
from cobra.core.metadata.keyvaluepair import ListOfKeyValue


class MetaData(collections.MutableMapping):
    """Class representation of the meta-data of an object.

    Meta-data consists of three components:
    - CVTerms, storing resource:identifier annotation information (this is
        exposed via the dict interface)
    - History, storing the object history
    - KeyValuePairs, a list of key-value pairs

    Parameters
    ----------
    cvterms : dict or CVTerms object
        The cvterm holds data for external resources
    history : dict, History
        The history is holding the data about the creator,
        created and modified dates.
    keyvalue_data : list
        Some key-value pairs which are not suitable to be
        represented anywhere else in the model.
    """

    def __init__(
        self,
        cvterms: Union[Dict, "CVTerms"] = None,
        history: Union[Dict, "History"] = None,
        keyvalue_data: List = None,
    ):
        self._cvterms = CVTerms()
        self.add_cvterms(cvterms)
        self.history = History.parse_history(history)
        self.key_value_data = ListOfKeyValue.parse_listofkeyvalue(keyvalue_data)

    def add_cvterm(self, cvterm, index: int = 0) -> None:
        """Adds a single CVTerm to CVList of given qualifier at given index."""
        self._cvterms.add_cvterm(cvterm=cvterm, index=index)

    def add_cvterms(self, cvterms: Union[Dict, "CVTerms"] = None) -> None:
        """Adds multiple CVTerm data."""
        self._cvterms.add_cvterms(cvterms=cvterms)

    @property
    def cvterms(self) -> "CVTerms":
        return self._cvterms

    @cvterms.setter
    def cvterms(self, value: Union[Dict, "CVTerms"]) -> None:
        self._cvterms = CVTerms()

        # synchronize the annotation dictionary
        self.cvterms._annotations = defaultdict(list)
        self.add_cvterms(value)

    @property
    def annotations(self) -> Dict:
        # annotation data in old format
        return self.cvterms.annotations

    def __getitem__(self, key: str) -> None:
        if key == "sbo" and len(self.annotations[key]) == 0:
            value = self._cvterms._annotations[key]
            del self._cvterms._annotations[key]
            return value
        return self.annotations[key]

    def __setitem__(self, key: str, value: List) -> None:
        self._cvterms.add_simple_annotations(dict({key: value}))

    def __delitem__(self, key: str) -> None:
        del self.annotations[key]

    def __iter__(self) -> Iterator:
        return iter(self.annotations)

    def __len__(self) -> int:
        return len(self.annotations)

    def __str__(self) -> str:
        return str(dict(self.annotations))

    def __repr__(self) -> str:
        return f"{dict(self.annotations)}"
