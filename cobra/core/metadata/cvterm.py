# -*- coding: utf-8 -*-

"""
Define the Controlled Vocabulary term class for refering to external
resources

"""
from __future__ import absolute_import

import re
import warnings
from collections import defaultdict
from collections.abc import MutableMapping, MutableSequence
from enum import Enum


class Qualifier(Enum):
    bqb_is = 1
    bqb_hasPart = 2
    bqb_isPartOf = 3
    bqb_isVersionOf = 4
    bqb_hasVersion = 5
    bqb_isHomologTo = 6
    bqb_isDescribedBy = 7
    bqb_isEncodedBy = 8
    bqb_encodes = 9
    bqb_occursIn = 10
    bqb_hasProperty = 11
    bqb_isPropertyOf = 12
    bqb_hasTaxon = 13
    bqb_unknown = 14
    bqm_is = 15
    bqm_isDescribedBy = 16
    bqm_isDerivedFrom = 17
    bqm_isInstanceOf = 18
    bqm_hasInstance = 19
    bqm_unknown = 20

URL_IDENTIFIERS_PATTERN = re.compile(
    r"^https?://identifiers.org/(.+?)[:/](.+)")


class CVTerm(object):
    """Representation of a single CVTerm."""
    def __init__(self, qualifier: Qualifier = Qualifier.bqb_is, resource: str = None):
        self.uri = resource
        if isinstance(qualifier, Qualifier):
            self.qualifier = qualifier
        else:
            raise TypeError("qualifier passed must be an enum Qualifier")

    def parse_provider_identifier(self):
        """Parses provider and term from given identifiers annotation uri.

        Parameters
        ----------
        uri : str
            uri (identifiers.org url)

        Returns
        -------
        (provider, identifier) if resolvable, None otherwise
        """
        match = URL_IDENTIFIERS_PATTERN.match(self.uri)
        if match:
            provider, identifier = match.group(1), match.group(2)
            if provider.isupper():
                identifier = "%s:%s" % (provider, identifier)
                provider = provider.lower()
        else:
            print("WARNING : %s does not conform to "
                  "'http(s)://identifiers.org/collection/id' or"
                  "'http(s)://identifiers.org/COLLECTION:id, so "
                  "is not added to annotation dictionary." % self.uri)
            return None

        return provider, identifier


# FIXME: this is probably not a dictionary
class CVTerms(MutableMapping):
    """
    Representation of all CVTerms of an object in their
    dependency structure.
    """

    # add property: annotations

    def __init__(self, data: dict = None):
        self._annotations = defaultdict(list)

        self._cvterms = defaultdict(list)
        if data is None:
            return
        elif isinstance(data, dict):
            for key, value in data.items():
                if key not in Qualifier.__members__:
                    raise TypeError("%s is not an enum Qualifier" % key)
                if isinstance(value, list):
                    self._cvterms[key] = CVList(value)
                else:
                    raise TypeError("The value passed must be of type list: "
                                    "{}".format(value))
        else:
            raise TypeError("Invalid format for CVTerms: '{}'".format(data))

    @staticmethod
    def parse_cvterms(data) -> 'CVTerms':
        """Tries to parse the CVterms."""
        if data is None or isinstance(data, dict):
            return CVTerms(data)
        elif isinstance(data, CVTerms):
            return data
        else:
            raise TypeError("Invalid format for CVTerms: '{}'".format(data))

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

        if index < len(self[qual]):
            self[qual][index]["resources"].append(cvterm.uri)
        elif index == len(self[qual]):
            self[qual].append({"resources":[cvterm.uri]})
        else:
            raise UnboundLocalError("The index is out of bound: {}".format(index))

    def add_cvterms(self, cvterms: 'CVTerms' = None):
        if cvterms is None:
            return
        elif isinstance(cvterms, dict) or isinstance(cvterms, CVTerm):
            parsed_cvterms = CVTerms.parse_cvterms(cvterms)
            for key, value in parsed_cvterms.items():
                # FIXME: this is probably broken now
                offset = len(self[key])
                for index in range(len(value)):
                    ex_res_list = value[index]
                    res_list = ex_res_list["resources"]  # FIXME: change to dot syntax
                    for uri in res_list:
                        cvterm = CVTerm(Qualifier[key], uri)
                        self.add_cvterm(cvterm, index+offset)
                    if "nested_data" in ex_res_list:
                        self[key][index+offset]["nested_data"] = ex_res_list["nested_data"]  # FIXME: change to dot syntax
        else:
            raise TypeError("The value passed must be of "
                            "type CVTerms: {}".format(cvterms))

    def __getitem__(self, key):
        return self._cvterms[key]

    def __setitem__(self, key, value):
        raise TypeError("The setitem method does not work for CVTerms. "
                        "Please use 'add_cvterm' method for adding cvterms.")

    def __delitem__(self, key):
        del self._cvterms[key]

    def __iter__(self):
        return iter(self._cvterms)

    def __len__(self):
        return len(self._cvterms)

    def __str__(self):
        return str(dict(self._cvterms))

    def __repr__(self):
        return '{}'.format(self._cvterms)


class CVList(MutableSequence):
    """
    Class representation of all sets of resources and their nested
    annotation corresponding to a given qualifier. It have similar
    structure like that of a list but has only restricted type of
    entries (of type ExternalResources) within it
    CVList : [
                 {
                    "resources" : [],
                    "nested_data" : CVTerm
                 },
                 {
                     "resources" : [],
                     "nested_data" : CVTerm
                 },
                 ...
              ]

    Parameters
    ----------
    cvlist : list
        a list containing entries confirming to ExternalResources structure

    """
    def __init__(self, data: list = None):

        self._sequence = list()
        if data is None:
            data = []
        elif not isinstance(data, list):
            raise TypeError("The data passed must be "
                            "inside a list: '{}'".format(data))

        for item in data:
            if isinstance(item, dict):
                self._sequence.append(ExternalResources(item))
            else:
                raise TypeError("All items inside CVList must be of type "
                                "dict: {}".format(item))

    def __len__(self):
        return len(self._sequence)

    def __delitem__(self, index):
        del self._sequence[index]

    def insert(self, index, value):
        if isinstance(value, ExternalResources):
            self._sequence.insert(index, value)
        elif isinstance(value, dict):
            self._sequence.insert(index, ExternalResources(value))
        else:
            raise TypeError("The passed format for setting external"
                            " resources is invalid.")

    def append(self, value):
        if isinstance(value, ExternalResources):
            self._sequence.append(value)
        elif isinstance(value, dict):
            self._sequence.append(ExternalResources(value))
        else:
            raise TypeError("The passed format for setting external"
                            " resources is invalid.")

    def __setitem__(self, index, value):
        if isinstance(value, ExternalResources):
            self._sequence[index] = value
        elif isinstance(value, dict):
            self._sequence[index] = ExternalResources(value)
        else:
            raise TypeError("The passed format for setting external"
                            " resources is invalid.")

    def __getitem__(self, index):
        return self._sequence[index]

    def __str__(self):
        return str(self._sequence)

    def __repr__(self):
        return '{}'.format(self._sequence)


class ExternalResources(MutableMapping):
    """
    Class representation of a single set of resources and its nested
    annotation. Its a special type of dict with restricted keys and
    values

    Parameters
    ----------
    data : dict
        A dictionary containing the resources and nested annotation
        {
            "resources" : [],
            "nested_data" : CVTerms
         }

    Allowed Keys
    ----------
    "resources" : string
        for accessing the mapped resources
    "nested_data" : string
        for accessing the nested annotation data

    """

    ANNOTATION_KEYS = ['resources', 'nested_data']

    def __init__(self, data=None):
        if data is None:
            data = {}
        self._mapping = dict()
        if not isinstance(data, dict):
            raise TypeError("The value passed must be of type dict.")
        for key, value in data.items():
            if key == 'resources':
                if not isinstance(value, list):
                    raise TypeError("Resources must be put in a list")
                self._mapping[key] = value
            elif key == 'nested_data':
                if isinstance(value, CVTerm):
                    self._mapping[key] = value
                elif isinstance(value, dict):
                    self._mapping[key] = CVTerms(value)
                else:
                    raise TypeError("The nested data structure does "
                                    "not have valid CVTerm format")
            elif key in Qualifier.__members__:
                self._mapping['nested_data'] = CVTerms({key: value})
            else:
                raise ValueError("Key '%s' is not allowed. Only "
                                 "allowed keys are 'resources', "
                                 "'nested_data'." % key)

    def __getitem__(self, key):
        if key not in self.ANNOTATION_KEYS:
            raise ValueError("Key %s is not allowed. Only allowed "
                             "keys are : 'resources', "
                             "'nested_data'" % key)
        return self._mapping[key]

    def __setitem__(self, key, value):
        """Restricting the keys and values that can be set.
           Only allowed keys are 'resources' and 'nested_data'
        """
        if key == 'resources':
            if not isinstance(value, list):
                raise TypeError("Resources must be put in a list")
            self._mapping[key] = value
        elif key == 'nested_data':
            if isinstance(value, CVTerm):
                self._mapping[key] = value
            elif isinstance(value, dict):
                self._mapping[key] = CVTerms(value)
            else:
                raise TypeError("The value passed has invalid format.")
        elif key in Qualifier.__members__:
            self._mapping['nested_data'] = CVTerms({key: value})
        else:
            raise ValueError("Key '%s' is not allowed. Only "
                             "allowed keys are 'resources', "
                             "'nested_data'." % key)

    def __delitem__(self, key):
        del self._mapping[key]

    def __iter__(self):
        return iter(self._mapping)

    def __len__(self):
        return len(self._mapping)

    def __str__(self):
        return str(self._mapping)

    def __repr__(self):
        return '{}'.format(self._mapping)

    def set_annotation(self, resource=None):
        if resource is None:
            return
        else:
            data = self._parse_annotation_info(resource)
            if data is None:
                return
            else:
                provider, identifier = data
            self.metadata["annotation"][provider].append((self.qualifier_key, identifier))
