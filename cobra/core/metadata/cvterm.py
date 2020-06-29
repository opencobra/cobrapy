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


# the supported qualifiers for cvterm
class Qualifier(Enum):
    bqb_is = 0
    bqb_hasPart = 1
    bqb_isPartOf = 2
    bqb_isVersionOf = 3
    bqb_hasVersion = 4
    bqb_isHomologTo = 5
    bqb_isDescribedBy = 6
    bqb_isEncodedBy = 7
    bqb_encodes = 8
    bqb_occursIn = 9
    bqb_hasProperty = 10
    bqb_isPropertyOf = 11
    bqb_hasTaxon = 12
    bqb_unknown = 13
    bqm_is = 14
    bqm_isDescribedBy = 15
    bqm_isDerivedFrom = 16
    bqm_isInstanceOf = 17
    bqm_hasInstance = 18
    bqm_unknown = 19

# the URL pattern to out parse provider and identifier
URL_IDENTIFIERS_PATTERN = re.compile(
    r"^https?://identifiers.org/(.+?)[:/](.+)")


class CVTerm(object):
    """Representation of a single CVTerm.

       Parameters
       ----------
       qualifier : Qualifier
            the qualifier relation of resource to the component
       resource : string
            a uri identifying external resource

        Attributes
        ----------
        qualifier : Qualifier
             the qualifier relation of resource to the component
        uri : string
             a uri identifying external resource
    """

    def __init__(self, qualifier: 'Qualifier' = Qualifier.bqb_is,
                 resource: 'str' = None):
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


class CVTerms(MutableMapping):
    """
    Representation of all CVTerms of an object in their
    dependency structure. It is like a dictionary where
    qualifier will be keys and CVList will be corresponding values

    Parameters
    ----------
    data : dict
        a dictionary mapping qualifier to its CVList/List

    This is how a CVTerm looks :
    {
        "bqb_is": [
            {
                "resources": [
                    "resource_uri",
                    ...
                ],
                "nested_data": CVTerms Object
            },
            ...
        ],
        ...
    }
    """

    def __init__(self, data: 'dict' = None):
        self._annotations = defaultdict(list)
        self._cvterms = defaultdict(CVList)
        if data is None:
            return
        elif isinstance(data, dict):
            for key, value in data.items():
                if key not in Qualifier.__members__:
                    raise TypeError("%s is not an enum Qualifier" % key)
                if isinstance(value, CVList):
                    self._cvterms[key] = value
                elif isinstance(value, list):
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
        """
        Adds a single CVTerm to CVTerms.

        Parameters
        ----------
        cvterm : CVTerm
            the cvterm to be added
        index : int
            the index where this cvterm should be added inside
            the CVList of corresponding qualifier
        """
        if isinstance(cvterm, CVTerm):
            qual = str(cvterm.qualifier)
            qual = qual[10:] if qual.startswith('Qualifier.') else qual
            data = cvterm.parse_provider_identifier()
            if data is not None:
                provider, identifier = data
                self._annotations[provider].append(identifier)
        else:
            raise TypeError("The CVTerm passed must be a CVTerm "
                            "object: {}".format(cvterm))

        if index < len(self[qual]):
            self[qual][index].resources.append(cvterm.uri)
        elif index == len(self[qual]):
            self[qual].append({"resources": [cvterm.uri]})
        else:
            raise UnboundLocalError("The index is out of bound:"
                                    " {}".format(index))

    def add_cvterms(self, cvterms: 'CVTerms' = None):
        """
        Adds multiple CVTerm to CVTerms.

        Parameters
        ----------
        cvterm : CVTerms or dict (to be converted in CVTerms)
            the cvterms to be added
        """
        if cvterms is None:
            return
        elif isinstance(cvterms, dict) or isinstance(cvterms, CVTerms):
            parsed_cvterms = CVTerms.parse_cvterms(cvterms)
            for key, value in parsed_cvterms.items():

                offset = len(self[key])
                for index in range(len(value)):
                    external_res = value[index]
                    res_list = external_res.resources
                    for uri in res_list:
                        cvterm = CVTerm(Qualifier[key], uri)
                        self.add_cvterm(cvterm, index+offset)
                    if external_res.nested_data is not None and \
                            len(external_res.nested_data) != 0:
                        self[key][index+offset].nested_data = \
                            external_res.nested_data
        else:
            raise TypeError("The value passed must be of "
                            "type CVTerms: {}".format(cvterms))

    def add_simple_annotations(self, data: None):
        """
        Adds cvterms via old annotation format. If no qualifier
        is linked to the identifier, default qualifier i.e "bqb_is"
        will be used.

        Parameters
        ----------
        data : dict
            the data in old annotation format
        """
        if data is None:
            data = {}

        # if annotation is in the form of list of list, modify the format
        if isinstance(data, list):
            dict_anno = defaultdict(list)
            for item in data:
                cvt = CVTerm(resource=item[1])
                data = cvt.parse_provider_identifier()
                if data is None:
                    continue
                else:
                    provider, identifier = data

                dict_anno[provider].append(identifier)
            data = dict_anno

        if not isinstance(data, dict):
            raise TypeError("The data passed must be of type "
                            "dict: {}".format(data))

        for key, value in data.items():

            # addition of "sbo" term
            if key == "sbo":
                if isinstance(value, str):
                    self._annotations[key] = list([value])
                elif isinstance(value, list):
                    self._annotations[key] = list(value)
                else:
                    raise TypeError("'sbo' terms must be wrapped"
                                    " inside a list: {}".format(value))
                continue

            # if single identifiers are put directly as string,
            # put them inside a list
            # for eg:
            # { "chebi": "CHEBI:17234"} -> { "chebi": ["CHEBI:17234"]}
            if isinstance(value, str) and key != 'sbo':
                data[key] = [value]
                value = [value]
            if not isinstance(value, list):
                raise TypeError("The value passed must be of type "
                                "list: {}".format(value))
            if not isinstance(key, str):
                raise TypeError("The key passed must be of type "
                                "string: {}".format(key))

            # reset the data of annotations corresponding to this key
            self._annotations[key] = []
            for identifier in value:
                cvterm = CVTerm()
                # if no qualifier is linked to identifier i.e annotation
                # of the form { "chebi": ["CHEBI:17234"]}
                if isinstance(identifier, str):
                    cvterm.uri = "https://identifiers.org/" + key + \
                                 "/" + identifier
                    cvterm.qualifier = Qualifier["bqb_is"]
                # if some qualifier is linked to the identifier i.e annotation
                # of the form { "chebi": ["bqb_is", "CHEBI:17234"]}
                elif isinstance(identifier, list):
                    cvterm.uri = "https://identifiers.org/" + key + "/" \
                                  + identifier[1]
                    cvterm.qualifier = Qualifier[identifier[0]]
                else:
                    raise TypeError("The identifier passed must be of \
                                     type string: {}".format(identifier))
                self.add_cvterm(cvterm, 0)

    def equals(self, cvterms):
        """
        Compare two CVTerms objects to find out whether they
        are same (have same data) or not
        """
        if len(self._cvterms) != len(cvterms):
            return False
        for key, value in cvterms.items():
            if key not in self._cvterms:
                return False
            if not value.equals(self._cvterms[key]):
                return False
        return True

    @property
    def annotations(self):
        return getattr(self, "_annotations", defaultdict(list))

    def __getitem__(self, key):
        if key not in Qualifier.__members__:
            raise TypeError("''%s' is not an valid enum Qualifier" % key)
        return self._cvterms[key]

    def __setitem__(self, key, value):
        """Make sure that key passed is of type string and value
           passed confirms to CVList type (CVList or list)
        """
        # setting the cvterm
        if key not in Qualifier.__members__:
            raise TypeError("%s is not an enum Qualifier" % key)
        if isinstance(value, list):
            self._cvterms[key] = CVList(value)
        elif isinstance(value, CVList):
            self._cvterms[key] = value
        else:
            raise TypeError("The value passed must be of type list"
                            " or CVList: {}".format(value))
        # setting the annotation
        for ex_res in value:
            for uri in ex_res.resources:
                cvterm = CVTerm(Qualifier[key], uri)
                data = cvterm.parse_provider_identifier()
                if data is not None:
                    provider, identifier = data
                    self._annotations[provider].append(identifier)

    def __delitem__(self, key):
        del self._cvterms[key]

    def __iter__(self):
        return iter(self._cvterms)

    def __len__(self):
        return len(self._cvterms)

    def __str__(self):
        return str(dict(self._cvterms))

    def __repr__(self):
        return '{}'.format(dict(self._cvterms))


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
    def __init__(self, data: 'list' = None):

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
                raise TypeError("All items inside CVList must be "
                                "of type dict: {}".format(item))

    def equals(self, cvlist):
        """
        Compare two CVList objects to find out whether
        they are same (have same data) or not
        """
        if len(self) != len(cvlist):
            return False
        num_ext_res = len(self)
        for index in range(num_ext_res):
            if not self[index].equals(cvlist[index]):
                return False
        return True

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
        return '{}'.format(list(self._sequence))


class ExternalResources(object):
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

    def __init__(self, data=None):
        if data is None:
            data = {}
        if not isinstance(data, dict):
            raise TypeError("The value passed must be of type dict.")
        for key, value in data.items():
            if key == 'resources':
                if not isinstance(data["resources"], list):
                    raise TypeError("Resources must be wrapped ""in a list:"
                                    " {}".format(data["resources"]))
                else:
                    self._resources = data["resources"]
            elif key == 'nested_data':
                if isinstance(value, CVTerms):
                    self._nested_data = value
                elif isinstance(value, dict):
                    self._nested_data = CVTerms(value)
                else:
                    raise TypeError("The nested data structure does "
                                    "not have valid CVTerm format: "
                                    "{}".format(value))
            elif key in Qualifier.__members__:
                self._nested_data = CVTerms({key: value})
            else:
                raise ValueError("Key '%s' is not allowed. Only "
                                 "allowed keys are 'resources', "
                                 "'nested_data'." % key)

    def equals(self, ext_res):
        """
        Compare two ExternalResources objects to find out whether
        they are same (have same data) or not
        """
        if self.resources != ext_res.resources:
            return False
        if self.nested_data is None and ext_res.nested_data is None:
            return True
        elif self.nested_data is None or ext_res.nested_data is None:
            return False
        elif not self.nested_data.equals(ext_res.nested_data):
            return False
        return True

    @property
    def resources(self):
        return getattr(self, "_resources", None)

    @resources.setter
    def resources(self, value):
        if not isinstance(value, list):
            raise TypeError("The resources must be wrapped "
                            "inside a list: {}".format(value))
        else:
            self._resources = value

    @property
    def nested_data(self):
        return getattr(self, "_nested_data", None)

    @nested_data.setter
    def nested_data(self, value):
        if isinstance(value, CVTerms):
            self._nested_data = value
        elif isinstance(value, dict):
            self._nested_data = CVTerms(value)
        else:
            raise TypeError("The nested data structure does "
                            "not have valid CVTerm format: {}".format(value))

    def __str__(self):
        if self.nested_data is None:
            return str({"resources": self.resources})
        else:
            return str({
                "resources": self.resources,
                "nested_data": self.nested_data
                })

    def __repr__(self):
        return self.__str__()
