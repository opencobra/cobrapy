# -*- coding: utf-8 -*-

"""
Define the Controlled Vocabulary term class for refering to external
resources
"""

from __future__ import absolute_import

import re
import warnings
from collections.abc import MutableMapping, MutableSequence


URL_IDENTIFIERS_PATTERN = re.compile(
    r"^https?://identifiers.org/(.+?)[:/](.+)")


class CVTerm(MutableMapping):
    """
    Class representation of Controlled Vocabulary term inside Annotation.
    It will look similar to a dictionary, but will have restrictions on the
    type of keys and values

    Parameters
    ----------
    cvterm : dict
        A dictionary type structure that maps the provider to its
        corresponding CVList

    Attributes
    ----------
        The providers are set as keys for accessing their CVLists

    """

    def __init__(self, metadata=None, cvterm=None):
        if cvterm is None:
            cvterm = {}
        self._mapping = dict()
        self.metadata = metadata
        if not isinstance(cvterm, dict):
            raise TypeError("The annotation data must be in a dict form")
        else:
            for key, value in cvterm.items():
                if not isinstance(key, str):
                    raise TypeError("the provider must be of type string")
                if isinstance(value, list):
                    self._mapping[key] = self.CVList(self.metadata, key, value)
                elif isinstance(value, self.CVList):
                    self._mapping[key] = value
                else:
                    raise TypeError("the value passed for key '%s' "
                                    "has invalid format" % key)

    def __getitem__(self, key):
        return self._mapping[key]

    def __setitem__(self, key, value):
        """Make sure that key passed is of type string and value
           passed confirms to CVList type (CVList or list)
        """
        if not isinstance(key, str):
            raise TypeError("The key passed must be a string")
        if isinstance(value, list):
            self._mapping[key] = self.CVList(self.metadata, key, value)
        elif isinstance(value, self.CVList):
            self._mapping[key] = value
        else:
            raise TypeError("The value passed does not confirm to CVList type")

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

    class CVList(MutableSequence):
        """
        Class representation of all sets of resources and their nested
        annotation corresponding to a given qualifier. It have similar
        structure like that of a list but has only restricted type of
        entries (of type ExternalResources) within it

        Parameters
        ----------
        cvlist : list
            a list containing entries confirming to ExternalResources structure

        """

        def __init__(self, metadata=None, qualifier_key=None, cvlist=None):
            if cvlist is None:
                cvlist = []
            if qualifier_key is None:
                self._qualifier_key = "is"
            elif not isinstance(qualifier_key, str):
                raise ("The qualifier key passed must be of type string")
            else:
                self._qualifier_key = qualifier_key
            self._sequence = list()
            self.metadata = metadata
            if not isinstance(cvlist, list):
                raise TypeError("The resources passed must be inside a list")
            for item in cvlist:
                if isinstance(item, CVTerm.ExternalResources):
                    self._sequence.append(item)
                elif isinstance(item, dict):
                    self._sequence.append(CVTerm.ExternalResources(self.metadata, self._qualifier_key, item))
                else:
                    raise TypeError("All items must confirm to "
                                    "ExternalResources structure")

        def __len__(self):
            return len(self._sequence)

        def __delitem__(self, index):
            del self._sequence[index]

        def insert(self, index, value):
            if isinstance(value, CVTerm.ExternalResources):
                self._sequence.insert(index, value)
            elif isinstance(value, dict):
                self._sequence.insert(index, CVTerm.ExternalResources(self.metadata, self._qualifier_key, value))
            else:
                raise TypeError("The passed format for setting external"
                                " resources is invalid.")

        def append(self, value):
            if isinstance(value, CVTerm.ExternalResources):
                self._sequence.append(value)
            elif isinstance(value, dict):
                self._sequence.append(CVTerm.ExternalResources(self.metadata, self._qualifier_key, value))
            else:
                raise TypeError("The passed format for setting external"
                                " resources is invalid.")

        def __setitem__(self, index, value):
            if isinstance(value, CVTerm.ExternalResources):
                self._sequence[index] = value
            elif isinstance(value, dict):
                self._sequence[index] = CVTerm.ExternalResources(self.metadata, self._qualifier_key, value)
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
                "nested_data" : CVTerm
             }

        Allowed Keys
        ----------
        "resources" : string
            for accessing the mapped resources
        "nested_data" : string
            for accessing the nested annotation data

        """

        ANNOTATION_KEYS = ['resources', 'nested_data', 'qualifier_type']
        QUALIFIER_RELATION = ['MODEL', 'BIOLOGICAL', 'UNKNOWN']

        def __init__(self, metadata=None, qualifier_key=None, data=None):
            if data is None:
                data = {}
            if qualifier_key is None:
                self.qualifier_key = "is"
            elif not isinstance(qualifier_key, str):
                raise TypeError("The qualifier key passed must be of type string")
            else:
                self.qualifier_key = qualifier_key
            self._mapping = dict()
            self.metadata = metadata
            if not isinstance(data, dict):
                raise TypeError("The value passed must be of type dict.")
            for key, value in data.items():
                if key not in self.ANNOTATION_KEYS:
                    raise ValueError("Key '%s' is not allowed. Only "
                                     "allowed keys are 'resources', "
                                     "'nested_data'." % key)
                if key == 'resources':
                    if not isinstance(value, list):
                        raise TypeError("Resources must be put in a list")
                    self._mapping[key] = value
                    for items in value:
                        self.set_annotation(items)
                elif key == 'nested_data':
                    if isinstance(value, CVTerm):
                        self._mapping[key] = value
                    elif isinstance(value, dict):
                        self._mapping[key] = CVTerm(value)
                    else:
                        raise TypeError("The nested data structure does "
                                        "not have valid CVTerm format")
                elif key == "qualifier_type":
                    if not isinstance(value, int):
                        raise TypeError("The value passed for qualifier type "
                                        "must be an integer")
                    if value == 0 or value == 1:
                        self._mapping[key] = self.QUALIFIER_RELATION[value]
                    else:
                        self._mapping[key] = self.QUALIFIER_RELATION[2]

        def __getitem__(self, key):
            if key not in self.ANNOTATION_KEYS:
                raise ValueError("Key %s is not allowed. Only allowed "
                                 "keys are : 'qualifier_type', 'resources', "
                                 "'nested_data'" % key)
            return self._mapping[key]

        def __setitem__(self, key, value):
            """Restricting the keys and values that can be set.
               Only allowed keys are 'resources' and 'nested_data'
            """
            if key not in self.ANNOTATION_KEYS:
                raise ValueError("Key %s is not allowed. Only allowed "
                                 "keys are : 'qualifier_type', 'resources', "
                                 "'nested_data'" % key)
            if key == 'resources':
                if not isinstance(value, list):
                    raise TypeError("Resources must be put in a list")
                self._mapping[key] = value
                for items in value:
                    set_annotation(items)
            elif key == 'nested_data':
                if isinstance(value, CVTerm):
                    self._mapping[key] = value
                elif isinstance(value, dict):
                    self._mapping[key] = CVTerm(value)
                else:
                    raise TypeError("The value passed has invalid format.")
            elif key == "qualifier_type":
                if not isinstance(value, int):
                    raise TypeError("The value passed for qualifier type "
                                    "must be an integer")
                if value == 0 or value == 1:
                    self._mapping[key] = self.QUALIFIER_RELATION[value]
                else:
                    self._mapping[key] = self.QUALIFIER_RELATION[2]

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

        def _parse_annotation_info(self, uri):
            """Parses provider and term from given identifiers annotation uri.

            Parameters
            ----------
            uri : str
                uri (identifiers.org url)

            Returns
            -------
            (provider, identifier) if resolvable, None otherwise
            """
            match = URL_IDENTIFIERS_PATTERN.match(uri)
            if match:
                provider, identifier = match.group(1), match.group(2)
                if provider.isupper():
                    identifier = "%s:%s" % (provider, identifier)
                    provider = provider.lower()
            else:
                print("WARNING : %s does not conform to "
                      "'http(s)://identifiers.org/collection/id' or"
                      "'http(s)://identifiers.org/COLLECTION:id, so "
                      "is not added to annotation dictionary." % uri)
                return None

            return provider, identifier
