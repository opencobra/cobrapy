"""Define the Controlled Vocabulary term class."""

import collections
import re
from collections import defaultdict
from enum import Enum
from typing import Dict, FrozenSet, Iterator, List, Optional, Tuple, Union
from warnings import warn

from .helper import URL_IDENTIFIERS_PATTERN, _parse_identifiers_uri


class Qualifier(Enum):
    """The possible qualifiers inside a CVTerm"""

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


class CVTerm:
    """Representation of a single CVTerm.

    Parameters
    ----------
    qualifier : Qualifier
         the qualifier relation of resource to the component
    resource : string
         a uri identifying external resource
    """

    def __init__(self, qualifier: Qualifier = Qualifier.bqb_is, resource: str = None):
        self.uri = resource
        if isinstance(qualifier, Qualifier):
            self.qualifier = qualifier
        elif isinstance(qualifier, str):
            if qualifier not in Qualifier.__members__:
                raise TypeError(f"{qualifier} is not a supported enum Qualifier")
            else:
                self.qualifier = Qualifier[qualifier]
        else:
            raise TypeError(f"{qualifier} is not a supported enum Qualifier")

    def parse_provider_identifier(self) -> Optional[Tuple]:
        """Parses provider and identifier term from given resource uri.

        Returns
        -------
        (provider, identifier) if resolvable, None otherwise
        """
        if self.uri is None:
            raise ValueError(f"'uri' set for this cvterm is None: {self}")
        else:
            return _parse_identifiers_uri(self.uri)


class CVTerms(collections.MutableMapping):
    """
    Representation of all CVTerms of an object in their
    dependency structure. It is like a dictionary where
    qualifier will be keys and CVList will be corresponding values.

    Parameters
    ----------
    data : dict
        a dictionary mapping qualifier to its CVList/List

    This is how a CVTerms looks :
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

    1. The only way to add annotation data via old format is by
       using the method "add_simple_annotations()".
    2. Single CVTerm data can be added by using "add_cvterm()".
    3. Multiple CVTerm data can be added by using "add_cvterms()".
    """

    def __init__(self, data: Dict = None):

        # storing data in new annotation format
        # as described above
        self._cvterms = defaultdict(CVList)

        if data is None:
            return
        elif isinstance(data, dict):
            for key, value in data.items():
                self.__setitem__(key, value)
        else:
            raise TypeError(f"Invalid format for CVTerms: '{data}'")

    @staticmethod
    def from_data(data: Union[Dict, "CVTerms"]) -> "CVTerms":
        """Parses a CVTerms object from given data"""
        if data is None or isinstance(data, dict):
            return CVTerms(data)
        elif isinstance(data, CVTerms):
            return data
        else:
            raise TypeError(f"Invalid format for CVTerms: '{data}'")

    def to_dict(self) -> dict:
        """Represent a CVTerms object in python dict.

        Returns:
        -------
        dict:
            a dict where each key has a list of all external resources in the original
            self._cvterms dictionary for that key

        Note: If this should work on Python 3.5 and lower, it needs to be an
        OrderedDict.

        """
        return {
            key: [ex_res.to_dict() for ex_res in value]
            for key, value in self._cvterms.items()
        }

    def add_cvterm(self, cvterm: CVTerm, index: int = 0) -> None:
        """
        Adds a single CVTerm to CVTerms at given index position
        corresponding to the passed qualifier.

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
            qual = qual[10:] if qual.startswith("Qualifier.") else qual
        else:
            raise TypeError(f"The CVTerm passed must be a CVTerm object: {cvterm}")

        if index < len(self[qual]):
            self[qual][index].resources.append(cvterm.uri)
        elif index == len(self[qual]):
            self[qual].append({"resources": [cvterm.uri]})
        else:
            raise UnboundLocalError(f"The index is out of bound: {index}")

    def add_cvterms(self, cvterms: Union[Dict, "CVTerms"] = None) -> None:
        """
        Adds multiple CVTerm to CVTerms.

        Parameters
        ----------
        cvterms : CVTerms or dict (to be added in CVTerms dict)
            the cvterms to be added
        """
        if cvterms is None:
            return

        if isinstance(cvterms, dict) or isinstance(cvterms, CVTerms):
            parsed_cvterms = CVTerms.from_data(cvterms)
            for key, value in parsed_cvterms.items():

                offset = len(self[key])
                for index in range(len(value)):
                    external_res = value[index]
                    res_list = external_res.resources
                    for uri in res_list:
                        cvterm = CVTerm(Qualifier[key], uri)
                        self.add_cvterm(cvterm, index + offset)
                    if (
                        external_res.nested_data is not None
                        and len(external_res.nested_data) != 0
                    ):
                        self[key][index + offset].nested_data = external_res.nested_data
        else:
            raise TypeError(f"The value passed must be of type CVTerms: {cvterms}")

    def add_simple_annotations(self, data: Union[Dict, List] = None) -> None:
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

        if not isinstance(data, (dict, list)):
            raise TypeError(f"The data passed must be of type dict or list: {data}")

        for key, value in data.items():

            # if single identifiers are put directly as string,
            # put them inside a list. For eg:
            # { "chebi": "CHEBI:17234"} -> { "chebi": ["CHEBI:17234"]}
            if isinstance(value, str) and key != "sbo":
                data[key] = [value]
                value = [value]
            if not isinstance(value, (list, str)):
                raise TypeError(
                    f"The value passed must be of type list or str: {value}"
                )

            # adding data one by one
            for identifier in value:
                qual = Qualifier["bqb_is"]
                # if no qualifier is linked to identifier i.e annotation
                # of the form { "chebi": "CHEBI:17234"}
                if isinstance(identifier, str):
                    uri = "https://identifiers.org/" + key + "/" + identifier
                # if some qualifier is linked to the identifier i.e annotation
                # of the form { "chebi": ["bqb_is", "CHEBI:17234"]}
                elif isinstance(identifier, list):
                    uri = "https://identifiers.org/" + key + "/" + identifier[1]
                    qual = Qualifier[identifier[0]]
                else:
                    raise TypeError(
                        f"The identifier passed must be of type string "
                        f"or list: {identifier}"
                    )
                self.add_cvterm(CVTerm(qualifier=qual, resource=uri), 0)

    @property
    def annotations(self) -> Dict:
        annotation_dict = dict()
        resources = self.resources
        for res in resources:
            if re.match(URL_IDENTIFIERS_PATTERN, res):
                provider, identifier = _parse_identifiers_uri(res)
                if provider in annotation_dict.keys():
                    annotation_dict[provider].append(identifier)
                else:
                    annotation_dict[provider] = [identifier]
        return {k: sorted(annotation_dict[k]) for k in sorted(annotation_dict.keys())}

    @property
    def resources(self) -> FrozenSet:
        """Get all resources.

        Returns:
        -------
        FrozenSet:
            a set of all external resources in the original self._cvterms dictionary
            including external resources of nested data
        """
        resources = set()
        for value in self._cvterms.values():
            for ex_res in value:
                resources.update(ex_res.resources)
                if ex_res.nested_data:
                    resources.update(ex_res.nested_data.resources)
        return frozenset(resources)

    def __getitem__(self, key: str) -> "CVList":
        """Get CVList corresponding to a qualifier"""
        if key not in Qualifier.__members__:
            raise TypeError(f"'{key}' is not an valid Qualifier.")
        return self._cvterms[key]

    def __setitem__(self, key: str, value: Union[List, "CVList"]) -> None:
        """Adds data to the CVTerm dict. 'key' has to be a Qualifier
        name string like "bqb_is" or "bqm_is" and value has to be a
        List/CVList object corresponding to that qualifier.
        """
        # setting the cvterm
        if key not in Qualifier.__members__:
            raise TypeError(f"{key} is not an enum Qualifier")
        if isinstance(value, list):
            self._cvterms[key] = CVList(value)
        elif isinstance(value, CVList):
            self._cvterms[key] = value
        else:
            raise TypeError(f"The value passed must be of type list or CVList: {value}")

    def __eq__(self, other: "CVTerms") -> bool:
        """Compare two CVTerms objects to find out whether they
        are same (have same data) or not
        """
        if len(self._cvterms) != len(other):
            return False
        for key, value in other.items():
            if key not in self._cvterms:
                return False
            if value != self._cvterms[key]:
                return False
        return True

    def __delitem__(self, key: str) -> None:
        del self._cvterms[key]

    def __iter__(self) -> Iterator:
        return iter(self._cvterms)

    def __len__(self) -> int:
        return len(self._cvterms)

    def __str__(self) -> str:
        return str(dict(self._cvterms))

    def __repr__(self) -> str:
        return self.__str__()


class CVList(collections.UserList):
    """
    Class representation of all sets of resources and their nested
    annotation corresponding to a given qualifier. It is a list but is restricted to
    entries of type ExternalResources within it. This list will
    have more than one entry inside it in the case of alternative
    set of resources linked to a given qualifier (such as "BQB_IS").

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
    data : list
        a list containing entries confirming to ExternalResources structure

    """

    def __init__(self, data: List[Union[Dict, "ExternalResources"]] = None):
        super().__init__(data)
        self.data = [
            ExternalResources(d_item)
            for d_item in self.data
            if not isinstance(d_item, ExternalResources)
        ]

    def insert(self, index: int, value: Union[Dict, "ExternalResources"]) -> None:
        """Insert a ExternalResource object at given index."""
        self.data.insert(index, _check_External_Resources_type(value))

    def append(self, value: Union[Dict, "ExternalResources"]) -> None:
        """Append a ExternalResource object to this list."""
        self.data.append(_check_External_Resources_type(value))

    def __setitem__(self, index: int, value: Union[Dict, "ExternalResources"]) -> None:
        self.data[index] = _check_External_Resources_type(value)


class ExternalResources:
    """
    Class representation of a single set of resources and its nested
    annotation. It is a special type of dict with restricted keys and
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
    "resources" : list
        for accessing the mapped resources
    "nested_data" : CVTerms
        for accessing the nested annotation data
    Qualifier: enum
        Members of the Qualifier enum defined above
    """

    def __init__(self, data: Dict = None):
        if data is None:
            data = {}
        if not isinstance(data, dict):
            raise TypeError("The value passed must be of type dict.")
        self.resources = data["resources"] if "resources" in data else {}
        self.nested_data = data["nested_data"] if "nested_data" in data else {}
        for key, value in data.items():
            if key == "resources":
                continue
            elif key == "nested_data":
                continue
            elif key in Qualifier.__members__:
                self._nested_data = CVTerms({key: value})
            else:
                raise ValueError(
                    f"Key '{key}' is not allowed. Only "
                    f"allowed keys are 'resources', "
                    f"'nested_data'."
                )

    @property
    def resources(self) -> List:
        return self._resources

    @resources.setter
    def resources(self, value: List) -> None:
        if value is None:
            self._resources = None
        elif not isinstance(value, list):
            raise TypeError(f"The resources must be wrapped inside a list: {value}")
        else:
            self._resources = value

    @property
    def nested_data(self) -> CVTerms:
        return self._nested_data

    @nested_data.setter
    def nested_data(self, value: Union[Dict, CVTerms]):
        if value is None:
            self._nested_data = None
        elif isinstance(value, CVTerms):
            self._nested_data = value
        elif isinstance(value, dict):
            self._nested_data = CVTerms(value)
        else:
            raise TypeError(
                f"The nested data structure does not have valid CVTerm format: {value}"
            )

    def to_dict(self):
        """Represents a ExternalResource object as python dict"""
        ex_dic = {"resources": list(self._resources)}
        if self.nested_data is not None:
            ex_dic["nested_data"] = self.nested_data.to_dict()
        return ex_dic

    def __eq__(self, other: "ExternalResources") -> bool:
        """
        Compare two ExternalResources objects to find out whether
        they are same (have same data) or not
        """
        if self.resources != other.resources:
            return False
        if self.nested_data is None and other.nested_data is None:
            return True
        elif self.nested_data is None or other.nested_data is None:
            return False
        elif self.nested_data != other.nested_data:
            return False
        return True

    def __str__(self) -> str:
        if self.nested_data is None:
            return str({"resources": self.resources})
        else:
            return str({"resources": self.resources, "nested_data": self.nested_data})

    def __repr__(self) -> str:
        return self.__str__()


def _check_External_Resources_type(
    value: Union[Dict, "ExternalResources"]
) -> "ExternalResources":
    if isinstance(value, dict):
        return ExternalResources(value)
    if isinstance(value, ExternalResources):
        return value
    raise TypeError(
        f"The passed object {value }for setting external resources has "
        f"invalid type: {type(value)}. It needs to be ExternalResouces or dict."
    )
