import collections
from typing import Dict, Iterator, List, Union

from cobra.core.metadata.cvterm import CVTerms
from cobra.core.metadata.history import History
from cobra.core.metadata.keyvaluepairs import KeyValuePairs


class MetaData(collections.MutableMapping):
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
        Some key-value pairs which are not suitable to be
        represented anywhere else in the model.
    """

    def __init__(
        self,
        cvterms: Union[Dict, CVTerms] = None,
        history: Union[Dict, History] = None,
        keyvaluepairs: List = None,
    ):
        self._cvterms = CVTerms.from_data(cvterms)
        self._history = History.from_data(history)
        self._keyvaluepairs = KeyValuePairs.from_data(keyvaluepairs)

    @property
    def annotations(self) -> Dict:
        """Backwards compatible annotations."""
        return self.cvterms.annotations

    def __setitem__(self, key: str, value: List) -> None:
        self._cvterms.add_simple_annotations(dict({key: value}))

    def __getitem__(self, key: str) -> None:
        '''
        if key == "sbo" and len(self.annotations[key]) == 0:
            # FIXME: what is this doing?
            value = self._cvterms._annotations[key]
            del self._cvterms._annotations[key]
            return value
        '''
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

    def add_cvterm(self, cvterm, index: int = 0) -> None:
        """Adds a single CVTerm to CVList of given qualifier at given index."""
        self.cvterms.add_cvterm(cvterm=cvterm, index=index)

    def add_cvterms(self, cvterms: Union[Dict, CVTerms] = None) -> None:
        """Adds multiple CVTerm data."""
        self.cvterms.add_cvterms(cvterms=cvterms)

    @property
    def history(self) -> History:
        return self._history

    @history.setter
    def history(self, history: Union[Dict, History]) -> None:
        self._history = History.from_data(history)
