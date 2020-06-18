# -*- coding: utf-8 -*-

from __future__ import absolute_import

from collections import defaultdict
from collections.abc import MutableMapping

from cobra.core.metadata.cvterm import CVTerms, CVTerm, Qualifier
from cobra.core.metadata.history import History
from cobra.core.metadata.keyvaluepair import ListOfKeyValue, KeyValueDict


class MetaData(MutableMapping):
    """Class representation of the meta-data of an object.

    Meta-data consists of three components:
    - CVTerms, storing resource:identifier annotation information (this is
        exposed via the dict interface)
    - History, storing the object history
    - KeyValuePairs, a dictionary of optional information

    Parameters
    ----------
    cvterms : dict or CVTerms object
        The cvterm holds data for external resources
    history : dict, History
        The history is holding the data about the creator,
        created and modified dates.
    keyValueDict : dict or KeyValuePair
        Some key-value pairs which are not suitable to be
        represented anywhere else in the model.
    """
    def __init__(self, cvterms: CVTerms = None, history: History = None,
                 keyValueDict: KeyValueDict = None):
        # internal dictionary of annotations for backwards compatibility
        # for resources a list of identifiers is stored

        self._cvterms = CVTerms()
        self.add_cvterms(cvterms)  # FIXME: unnecessary -> CVTerms(cvterms)
        self.history = History.parse_history(history)
        self.keyValueDict = KeyValueDict.parse_keyValueDict(keyValueDict)

    def add_cvterm(self, cvterm, index: int = 0):
        self._cvterms.add_cvterm(cvterm=cvterm, index=index)

    def add_cvterms(self, cvterms: CVTerms = None):
        self._cvterms.add_cvterms(cvterms=cvterms)

    @property
    def cvterms(self):
        return self._cvterms

    @cvterms.setter
    def cvterms(self, value):
        self._cvterms = CVTerms()
        
        # synchronize the annotation dictionary
        self.cvterms._annotations = defaultdict(list)
        self.add_cvterms(value)

    @property
    def annotations(self):
        return self.cvterms._annotations

    def __getitem__(self, key):
        return self.annotations[key]

    def __setitem__(self, key, value):
        self.annotations[key].append(value)

    def __delitem__(self, key):
        del self.annotations[key]

    def __iter__(self):
        return iter(self.annotations)

    def __len__(self):
        return len(self.annotations)

    def __str__(self):
        return str(dict(self.annotations))

    def __repr__(self):
        return '{}'.format(dict(self.annotations))
