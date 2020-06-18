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
        self._annotations = defaultdict(list)

        self._cvterms = CVTerms()
        self.add_cvterms(cvterms)
        self.history = History.parse_history(history)
        self.keyValueDict = KeyValueDict.parse_keyValueDict(keyValueDict)

    def add_cvterm(self, cvterm, index):
        if isinstance(cvterm, CVTerm):
            qual = str(cvterm.qualifier)
            qual = qual[10:] if qual.startswith('Qualifier.') else qual
            data = cvterm.parse_provider_identifier()
            if data is not None:
                provider, identifier = data
                self._annotations[provider].append(identifier)
        else:
            raise TypeError("The CVTerm passed must be a CVTerm object: {}".format(cvterm))

        if index < len(self._cvterms[qual]):
            self._cvterms[qual][index]["resources"].append(cvterm.uri)
        elif index == len(self._cvterms[qual]):
            self._cvterms[qual].append({"resources":[cvterm.uri]})
        else:
            raise UnboundLocalError("The index is out of bound: {}".format(index))

    def add_cvterms(self, cvterms: CVTerms = None):
        if cvterms is None:
            return
        elif isinstance(cvterms, dict) or isinstance(cvterms, CVTerm):
            parsed_cvterms = CVTerms.parse_cvterms(cvterms)
            for key, value in parsed_cvterms.items():
                offset = len(self.cvterms[key])
                for index in range(len(value)):
                    ex_res_list = value[index]
                    res_list = ex_res_list["resources"]
                    for uri in res_list:
                        cvterm = CVTerm(Qualifier[key], uri)
                        self.add_cvterm(cvterm, index+offset)
                    if "nested_data" in ex_res_list:
                        self._cvterms[key][index+offset]["nested_data"] = ex_res_list["nested_data"]
        else:
            raise TypeError("The value passed must be of "
                            "type CVTerms: {}".format(cvterms))

    @property
    def cvterms(self):
        return self._cvterms

    @cvterms.setter
    def cvterms(self, value):
        self._cvterms = CVTerms()
        self._annotations = defaultdict(list)
        self.add_cvterms(value)

    def __getitem__(self, key):
        return self._annotations[key]

    def __setitem__(self, key, value):
        self._annotations[key].append(value)

    def __delitem__(self, key):
        del self._annotations[key]

    def __iter__(self):
        return iter(self._annotations)

    def __len__(self):
        return len(self._annotations)

    def __str__(self):
        return str(dict(self._annotations))

    def __repr__(self):
        return '{}'.format(dict(self._annotations))
