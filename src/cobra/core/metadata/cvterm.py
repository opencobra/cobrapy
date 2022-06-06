"""Define the Controlled Vocabulary term class."""

import re
from collections import OrderedDict, UserList
from enum import Enum
from typing import Dict, FrozenSet, Iterator, List, Optional, Tuple, Union, Iterable, \
    Pattern, Callable, Any
import logging

from .helper import URL_IDENTIFIERS_PATTERN, parse_identifiers_uri

logger = logging.getLogger(__name__)


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
    """
    Representation of one CVTerm of an object in their
    dependency structure. It is a ????


    named tuple where the first named field is a qualifier and the second named field
    is ExternalResources.

    UserList keeps the data in the property data, which is a real list.


    Parameters
    ----------
    initialdata : dict
        a dictionary mapping qualifier to its CVList/List

    This is how a CVTerm looks :
    CVTerm.qualifier = "bqb_is"
    CVTerm.ex_res =
        {"resources": [
                    "resource_uri",
                    ...
                ],
                "nested_data": CVTerms Object
        }


    The internal dictionary is an ExternalResources

    1. The only way to add annotation data via old format is by
       using the method "add_simple_annotations()".
    2. Single CVTerm data can be added by using "add_cvterm()".
    3. Multiple CVTerm data can be added by using "add_cvterms()".
    """

    def __init__(
        self,
        ex_res: "ExternalResources" = None,
        qualifier: Union[Qualifier, str] = Qualifier.bqb_is,
    ):
        self._ex_res = None
        self._qualifier = None
        self._ex_res = self.check_ex_res_type(ex_res)
        self._qualifier = self.check_qualifier_type(qualifier)

    @property
    def qualifier(self) -> Qualifier:
        return self._qualifier

    @qualifier.setter
    def qualifier(self, qualifier: Union[str, int, Qualifier]) -> None:
        self._qualifier = self.check_qualifier_type(qualifier)

    @property
    def external_resources(self) -> "ExternalResources":
        return self._ex_res

    @external_resources.setter
    def external_resources(
        self, external_resources: Union[dict, "ExternalResources"]
    ) -> None:
        self._ex_res = self.check_ex_res_type(external_resources)

    @staticmethod
    def check_ex_res_type(
        ex_res: Union["ExternalResources", Dict]
    ) -> "ExternalResources":
        if ex_res is None:
            return ExternalResources()
        elif isinstance(ex_res, ExternalResources):
            return ex_res
        elif isinstance(ex_res, dict):
            return ExternalResources.from_dict(ex_res)
        else:
            raise TypeError(
                f"Allowed types for CVTerms ex_ress are ExternalResources "
                f"or dict, not {type(ex_res)}: {ex_res}"
            )

    @staticmethod
    def check_qualifier_type(qual: Union[int, str, Qualifier]) -> Qualifier:
        if (
            isinstance(qual, str)
            and qual not in Qualifier.__members__
            or (isinstance(qual, int) and not 0 <= qual < Qualifier.__len__())
        ):
            raise TypeError(f"{qual} is not a supported enum Qualifier")
        elif isinstance(qual, Qualifier):
            return qual
        elif isinstance(qual, str):
            return Qualifier[qual]
        elif isinstance(qual, int):
            return Qualifier(qual)
        else:
            raise TypeError(
                f"Allowed types for CVTerm qualifier must be Qualifier,"
                f"str member of the Qualifier enum "
                f"or an int between 0 and the Qualifier enum length "
                f"{type(qual)}, {qual}"
            )

    def to_ordered_dict(self) -> dict:
        """Represent a CVTerms object in python dict.

        Returns:
        -------
        dict:
            a dict where each key has is a qualifier and has a list of all external
            resources in the original self._cvterms dictionary for that key

        """
        return OrderedDict({self.qualifier: self.external_resources.to_dict()})

    @classmethod
    def from_dict(cls, data_dict: Dict):
        return cls(ex_res=data_dict.get("external_resources", None),
                   qualifier=data_dict.get("qualifier", Qualifier['bqb_is']))

    def __eq__(self, other: "CVTerm") -> bool:
        if not isinstance(other, CVTerm):
            return False
        if self.qualifier != other.qualifier:
            return False
        if self.external_resources != other.external_resources:
            return False
        return True


class CVTerms(UserList):
    """
    Representation of all CVTerms of an object in their
    dependency structure. It is list that contains qualifiers and external resouces
    for each CVTerm.

    Parameters
    ----------
    data : list
        a list containing qualifier and external resources in CVTerm forma

    This is how a CVTerms looks :
    [
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
    ]

    1. The only way to add annotation data via old format is by
       using the method "add_simple_annotations()".
    2. Single CVTerm data can be added by using "add_cvterm()".
    3. Multiple CVTerm data can be added by using "add_cvterms()".
    """

    def __init__(self, data=None):

        # storing data in new annotation format
        # as described above
        if data is None:
            data = []
        data = [self._check_CVTerm2(datum) for datum in data]
        super().__init__(data)

    @staticmethod
    def _check_CVTerm2(cvterm: Union["CVTerm", Dict]) -> Optional["CVTerm"]:
        if cvterm is None:
            return
        if isinstance(cvterm, CVTerm):
            return cvterm
        elif isinstance(cvterm, dict):
            return CVTerm.from_dict(cvterm)
        else:
            raise TypeError(
                f"Allowed types for CVTerms ex_ress are CVTerm "
                f"or dict, not {type(cvterm)}: {cvterm}"
            )

    @staticmethod
    def from_data(data: Union[List, Dict, "CVTerm"]) -> "CVTerms":
        """Parses a CVTerms object from given data"""
        if data is None:
            return CVTerms()
        if isinstance(data, dict):
            # TODO - need to check dict.py
            return CVTerms.from_dict(data)
        elif isinstance(data, list):
            return CVTerms(data)
        elif isinstance(data, CVTerm):
            return CVTerms([data])
        else:
            raise TypeError(f"Invalid format for CVTerms: '{data}'")

    def to_dict(self) -> dict:
        """Represent a CVTerms object in python dict.

        Returns:
        -------
        dict:
            a dict where each key has is a qualifier and has a list of all external
            resources in the original self.data list that match that qualifier

        """
        qualifier_set = set()
        for cvterm in self.data:
            qualifier_set.add(cvterm.qualifier.name)
        empty_lists = [[] for _ in qualifier_set]
        qualifier_dict = dict(zip(qualifier_set, empty_lists))
        for cvterm in self.data:
            qualifier_dict[cvterm.qualifier.name].append(cvterm.external_resources.to_dict())
        return OrderedDict(qualifier_dict)

    @classmethod
    def from_dict(cls, data_dict):
        cvterms_list = []
        for qual, ex_resources in data_dict.items():
            cvterms = [CVTerm(ex_res=ExternalResources.from_dict(ex_res),
                              qualifier=Qualifier[qual])
                       for ex_res in ex_resources]
            cvterms_list.extend(cvterms)
        return cls(cvterms_list)

    def add_cvterms(self, cvterms: Iterable[Union[Dict, "CVTerm"]]) -> None:
        """
        Adds multiple CVTerm to CVTerms.

        Parameters
        ----------
        cvterms : CVTerms list of CVTerm or dict (to be added in CVTerms dict)
            the cvterms to be added
        """
        self.extend(cvterms)

    def add_simple_annotations(self, data: Union[Dict, List] = None) -> None:
        """
        Adds cvterms via old annotation format. If no qualifier
        is linked to the identifier, default qualifier i.e "bqb_is"
        will be used.
        This function will add identifiers.org to as the URI.
        If the annotation does not match the identifier format, you should XXXXX?????

        Parameters
        ----------
        data : dict
            the data in old annotation format
        """
        if data is None:
            data = {}

        if not isinstance(data, (dict, list)):
            raise TypeError(f"The data passed must be of type dict or list: {data}")

        cvterm_list = []
        for key, value in data.items():

            # if single identifiers are put directly as string,
            # put them inside a list. For eg:
            # { "chebi": "CHEBI:17234"} -> { "chebi": ["CHEBI:17234"]}
            if isinstance(value, str):
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
                cvterm_list.append(CVTerm(ex_res=ExternalResources([uri]), qualifier=qual))
        self.add_cvterms(cvterm_list)

    @property
    def annotations(self) -> Dict:
        annotation_dict = {}
        resources = self.resources
        for res in resources:
            if re.match(URL_IDENTIFIERS_PATTERN, res):
                provider, identifier = parse_identifiers_uri(res)
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
        for datum in self.data:
            if datum.external_resources.resources:
                resources.update(datum.external_resources.resources)
            if datum.external_resources.nested_data:
                resources.update(datum.external_resources.nested_data.resources)
        return frozenset(resources)

    @property
    def qualifiers(self) -> FrozenSet:
        qualifier_set = set()
        for datum in self.data:
            qualifier_set.add(datum.qualifier)
        return frozenset(qualifier_set)

    def get_by_qualifier(self, qual_to_find: Union[str, Qualifier, int]) -> "CVTerms":
        qual_to_find = CVTerm.check_qualifier_type(qual_to_find)
        cvterms_to_return = [cvterm for cvterm in self.data
                             if cvterm.qualifier == qual_to_find]
        return self.__class__(cvterms_to_return)

    def __setitem__(self, key: int, value: CVTerm) -> None:
        UserList.__setitem__(self, key, self._check_CVTerm2(value))

    def append(self, item: CVTerm) -> None:
        """Append CVTerm to end."""
        self._check_CVTerm2(item)
        UserList.append(self, item)

    def extend(self, iterable: Union["CVTerms", Iterable[Union[CVTerm, Dict]]]) -> None:
        """Extend data list by appending elements from the iterable.

        Parameters
        ----------
        iterable : Iterable
        """
        if isinstance(iterable, CVTerms):
            self.data.extend(iterable.data)
        elif isinstance(iterable, Iterable):
            self.data.extend([self._check_CVTerm2(i) for i in iterable])

    def __eq__(self, other: "CVTerms") -> bool:
        """Compare two CVTerms objects to find out whether they
        are same (have same data) or not
        """
        if not isinstance(other, CVTerms):
            return False
        if len(self.data) != len(other.data):
            return False
        for self_i, other_i in zip(self.data, other.data):
            if self_i != other_i:
                return False
        return True

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
    """

    def __init__(self, resources: List = None, nested_data=None):
        self._resources = None
        self._nested_data = None
        self.resources = resources
        if resources:
            self.nested_data = nested_data

    @property
    def resources(self) -> List:
        return self._resources

    @resources.setter
    def resources(self, value: List) -> None:
        if value is None or isinstance(value, list):
            self._resources = value
        elif not isinstance(value, list):
            raise TypeError(f"The resources must be wrapped inside a list: {value}")

    @property
    def nested_data(self) -> CVTerms:
        return self._nested_data

    @nested_data.setter
    def nested_data(self, value: Union[Dict, CVTerm, CVTerms]):
        if value is None:
            self._nested_data = CVTerms()
        elif isinstance(value, CVTerm):
            self._nested_data = CVTerms([value])
        elif isinstance(value, CVTerms):
            self._nested_data = value
        elif isinstance(value, dict):
            self._nested_data = CVTerms.from_dict(value)
        else:
            raise TypeError(
                f"The nested data structure does not have valid CVTerm format: {value}"
            )

    def to_dict(self):
        """Represents a ExternalResource object as python dict"""
        ex_dic = {"resources": list(self._resources)}
        if self.nested_data is not None and len(self.nested_data):
            ex_dic["nested_data"] = self.nested_data.to_dict()
        return ex_dic

    @classmethod
    def from_dict(cls, input_data: Dict) -> "ExternalResources":
        """

        Parameters
        ----------
        input_data: dict

        Returns
        -------
        ExternalResources2: ExternalResources

        """
        ex_res = cls(
            resources=input_data.get("resources", None),
            nested_data=input_data.get("nested_data", None),
        )
        return ex_res

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
        return str(self.to_dict())

    def __repr__(self) -> str:
        return self.__str__()

    def _repr_html_(self):
        return f"""
                <p>
                    <strong>Resources</strong>{"<p>".join(self.resources)}<p>
                    <strong>Nested Data</strong>{self.nested_data}<p>
                    <strong>Memory address</strong>{id(self):#x}
                <p>"""
