"""Define the Controlled Vocabulary term class."""

from collections import UserList
from collections.abc import Iterable as ABCIterable
from enum import Enum
from typing import (
    Dict,
    Iterable,
    List,
    Optional,
    Union,
)
import re

from .helper import URL_IDENTIFIERS_PATTERN, parse_identifiers_uri


class Qualifier(Enum):
    """The possible qualifiers inside a CVTerm.

    The qualifiers and their detailed description are present in:
    https://co.mbine.org/author/biomodels.net-qualifiers/

    Qualifiers are divided into two groups
    Biological (bqb)    These kinds of qualifiers define the relationship between a
                        biological object represented by a model element and its
                        annotation.
    Modelling (bqm)     These kinds of qualifiers define the relationship between a
                        modelling object and its annotation.
    """

    def __init__(self, value):
        """Initialize Qualifier enums by creating a lookup dictionary."""
        self.__class__._map = getattr(self.__class__, "_map", {}) | {value: self}

    Biological_is = "bqb_is"
    Biological_hasPart = "bqb_hasPart"
    Biological_isPartOf = "bqb_isPartOf"
    Biological_isVersionOf = "bqb_isVersionOf"
    Biological_hasVersion = "bqb_hasVersion"
    Biological_isHomologTo = "bqb_isHomologTo"
    Biological_isDescribedBy = "bqb_isDescribedBy"
    Biological_isEncodedBy = "bqb_isEncodedBy"
    Biological_encodes = "bqb_encodes"
    Biological_occursIn = "bqb_occursIn"
    Biological_hasProperty = "bqb_hasProperty"
    Biological_isPropertyOf = "bqb_isPropertyOf"
    Biological_hasTaxon = "bqb_hasTaxon"
    Biological_unknown = "bqb_unknown"

    Modelling_is = "bqm_is"
    Modelling_isDescribedBy = "bqm_isDescribedBy"
    Modelling_isDerivedFrom = "bqm_isDerivedFrom"
    Modelling_isInstanceOf = "bqm_isInstanceOf"
    Modelling_hasInstance = "bqm_hasInstance"
    Modelling_unknown = "bqm_unknown"


class Identifier:
    def __init__(self, uri: str) -> None:
        self._uri = None
        self._namespace = None
        self._identifier = None

        self.uri = uri

    @property
    def uri(self) -> Optional[str]:
        return self._uri

    @uri.setter
    def uri(self, value: str) -> None:
        if re.match(URL_IDENTIFIERS_PATTERN, value):
            identifier_match = parse_identifiers_uri(value)

            if identifier_match is None:
                raise ValueError(f"The provided URI is not valid: {value}")
            namespace, identifier = identifier_match
            self._namespace = namespace
            self._identifier = identifier

    @property
    def namespace(self) -> Optional[str]:
        return self._namespace

    @property
    def identifier(self) -> Optional[str]:
        return self._identifier


class StandardizedAnnotation:
    """CVTerm class, representing controlled vocabulary.

    Controlled Vocabulary (CVTerm) can be defined as a curated and controlled
    relationship, described by a Qualifier (see above) - the relationship between an
    object and annotation must be part of the Qualifier class. These relationships
    are based in biochemical or biological relationships. The qualifiers/relationships
    are divided into bqbiol/bqb (biological qualification) and bqmodel/bqm (model
    qualifications). See two examples:
    "bqb_is" The biological entity represented by the SBML component is the subject
    of the referenced resource. This could serve to link a reaction to its counterpart
    in (e.g.) the ChEBI or Reactome databases.
    "bqm_is" The modeling object encoded by the SBML component is the subject of
    the referenced resource. This might be used, e.g., to link the model
    to an entry in a model database.
    See `Biomodels Qualifiers
    <https://co.mbine.org/author/biomodels.net-qualifiers/>`_
    For a definition of all qualifiers, see `SBML Level 3, Version 2 Core, p 104
    <https://identifiers.org/combine.specifications:sbml.level-3.version-2.core.release-2>`_

    The annotation will have one or more URI, which are encapsulated in
    ExternalResources class (see below).

    Each CVTerm has only ONE qualifier, and can have many URI in the resources. If
    you need to use another qualifier, it can be nested data (if relevant), or it
    should be in another CVTerm.
    If an object has multiple CVTerms, they are placed in a CVTermList (see below).

    This is how a CVTerm object looks :
    CVTerm.qualifier = "bqb_is"
    CVTerm.ex_res =
        {"resources": [
                    "resource_uri",
                    ...
                ],
                "nested_data":CVTermList Object
        }

    Examples of how CVTerms can be used

    Model examples (Each of these is a separate CVTerm)
    qualifier=bqm_is
    resources=["https://identifiers.org/biomodels.db/BIOMD0000000003"]
            A model identifier
    qualifier=bqm_isDescribed_by
    resources=["https://identifiers.org/pubmed/1833774"]
            A published article detailing the model
    qualifier=bqm_isVersionOf
    resources=["https://identifiers.org/wikipathways/WP179",
                "https://identifiers.org/reactome/REACT_152"/]
            Two links to what this model is a version of (in this case, cell cycle).

    Reaction examples
    qualifier=bqb_is
    resources=["https://identifiers.org/reactome/REACT_6327"/]
        A link to a reaction database that details reactions.
    qualifier=bqb_hasPart
    resources=["http://identifiers.org/uniprot/P04551",
                http://identifiers.org/uniprot/P10815"]
             resources.nested_date = {
                    qualifier=bqb_isDescribedby
                    resources=["https://identifiers.org/pubmed/1111111"]
        Two proteins that form part of the same complex. The nested data links to an
        article describing the formation of the complex.
        It is nested data because it is relevant to the hasPart CVTerm, but uses a
        different qualifier.
    """

    # TODO: Update
    def __init__(
        self,
        identifiers: Optional[
            Union["Identifier", str, Iterable[Union["Identifier", str]]]
        ] = None,
        qualifier: Union[Qualifier, str] = Qualifier.Biological_is,
        annotations: Optional[Iterable["StandardizedAnnotation"]] = None,
    ):
        """Initialize a CVTerm.

        Parameters
        ----------
        identifiers: Identifier or str or list
            The identifiers (URI format).
        qualifier: Qualifier or str
            The qualifier for the relationship.
        annotations: list
            List of StandardizedAnnotation objects that should be nested in this object.
        """
        self._identifiers = self.check_identifier_type(identifiers)
        self._qualifier = self.check_qualifier_type(qualifier)
        self._annotations = self.check_annotation_type(annotations)

    @property
    def qualifier(self) -> Qualifier:
        """Get qualifier for CVTerm.

        Returns
        -------
        Qualifier
        """
        return self._qualifier

    @qualifier.setter
    def qualifier(self, qualifier: Union[str, Qualifier]) -> None:
        """Set Qualifier.

        Parameters
        ----------
        qualifier - str, int or Qualifier
            Is converted to the Qualifier class.

        See Also
        --------
        CVTerm.check_qualifier_type()
        """
        self._qualifier = self.check_qualifier_type(qualifier)

    @property
    def identifiers(self) -> List["Identifier"]:
        """Get external identifiers.

        Returns
        -------
        ExternalResources
        """
        return self._identifiers

    @identifiers.setter
    def identifiers(self, identifiers: Iterable[Union[str, "Identifier"]]) -> None:
        """Set external resources.

        Parameters
        ----------
        identifiers - list of str or Identifier

        See Also
        --------
        CVTerm.check_identifier_type()
        """
        self._identifiers = self.check_identifier_type(identifiers)

    @property
    def annotations(self) -> List["StandardizedAnnotation"]:
        """Get the nested annotations.

        Returns
        -------
        List of annotations
        """
        return self._annotations

    @annotations.setter
    def annotations(self, annotations: Iterable["StandardizedAnnotation"]) -> None:
        """Set the nested annotations.

        Parameters
        ----------
        annotations - list of StandardizedAnnotation objects.

        See Also
        --------
        StandardizedAnnotation.check_annotation_type()
        """
        self._annotations = self.check_annotation_type(annotations)

    # @property
    # def resources(self) -> FrozenSet:
    #     """Get all resources.
    #
    #     Returns:
    #     -------
    #     FrozenSet:
    #         a set of all resources in the CVTerm as a set of strings
    #         including external resources of nested data as strings
    #     """
    #     return self.external_resources.resource_nested
    #
    @staticmethod
    def check_identifier_type(
        identifiers: Optional[
            Union["Identifier", str, Iterable[Union["Identifier", str]]]
        ],
    ) -> List["Identifier"]:
        """Check and parse input to ExternalResources.

        Parameters
        ----------
        ex_res: ExternalResources or dict or str, optional
            Input data to check if it is or can be transformed to ExternalResources
            class. String must start with http:// or https:// to be acceptable.
            Dictionary must match the format required by from_dict.
            If None is given, an empty ExternalResources is returned.
            No parsing of identifiers/URIs is done, perhaps in future versions.

        Returns
        -------
        ExternalResources

        Raises
        ------
        TypeError
            If given anything other than None, str, dict or ExternalResources.
            Will raise this error if given a string that does not start with http(s)://

        See Also
        --------
        ExternalResources.from_dict()
        """
        # TODO: Fix doc
        if identifiers is None:
            return []
        elif isinstance(identifiers, Identifier):
            return [identifiers]
        elif isinstance(identifiers, str):
            return [Identifier(identifiers)]
        elif isinstance(identifiers, ABCIterable):
            return [x for y in identifiers for x in __class__.check_identifier_type(y)]
        else:
            raise TypeError(
                f"Allowed types for identifiers are Identifier, str, or a list thereof,"
                f"not {type(identifiers)}: {identifiers}"
            )

    @staticmethod
    def check_qualifier_type(qual: Union[str, Qualifier]) -> Qualifier:
        """Check and parse input to Qualifier class.

        Parameters
        ----------
        qual: str or Qualifier, optional
            Input data to check if it is or can be transformed to Qualifier class.
            Strings must be a member of the Qualifier values.
            If None is given, an empty Qualifier is returned.

        Returns
        -------
        Qualifier

        Raises
        ------
        TypeError
            If given anything other than None, str, or Qualifier.
            Will raise this error if given a string that does not match the defined
            Qualifier members.
        """
        if isinstance(qual, str) and qual not in Qualifier._map:
            raise TypeError(f"{qual} is not a supported enum Qualifier")
        elif isinstance(qual, Qualifier):
            return qual
        elif isinstance(qual, str):
            return Qualifier._map[qual]
        else:
            raise TypeError(
                f"Allowed types for StandardAnnotation qualifiers are Qualifier or"
                f"str member of the Qualifier enum {type(qual)}, {qual}"
            )

    @staticmethod
    def check_annotation_type(
        ann: Optional[
            Union["StandardizedAnnotation", Iterable["StandardizedAnnotation"]]
        ],
    ) -> List["StandardizedAnnotation"]:
        if ann is None:
            return []
        elif isinstance(ann, StandardizedAnnotation):
            return [ann]
        elif isinstance(ann, ABCIterable):
            return [x for y in ann for x in __class__.check_annotation_type(y)]
        else:
            raise TypeError(
                f"Allowed types for nested StandardizedAnnotation annotations"
                f"are StandardizedAnnotation or a list of StandardizedAnnotation"
                f"objects, not {type(ann)}: ann"
            )

    # def to_dict(self) -> Dict:
    #     """Represent a CVTerm object in python dict.
    #
    #     Returns
    #     -------
    #     dict:
    #         A dict that has two keys
    #         "qualifier" - the qualifier as a string
    #         "external_resources" - the resources as a dictionary
    #
    #     See Also
    #     --------
    #     ExternalResources.to_dict()
    #
    #     """
    #     return {
    #         "qualifier": self.qualifier.value,
    #         "external_resources": self.external_resources.to_dict(),
    #     }
    #
    # @classmethod
    # def from_dict(cls, data_dict: Dict) -> "CVTerm":
    #     """Generate a CVTerm object based on a python dict.
    #
    #     Parameters
    #     ----------
    #     data_dict: dict
    #         A dict that has two keys
    #         "qualifier" - the qualifier as a string, optional. If not present, the
    #         qualifier is set to bqb_is.
    #         "external_resources" - the resources as a dictionary, optional
    #
    #     Returns
    #     -------
    #     CVTerm
    #
    #     See Also
    #     --------
    #     ExternalResources.to_dict()
    #
    #     """
    #     return cls(
    #         ex_res=data_dict.get("external_resources", None),
    #         qualifier=data_dict.get("qualifier", Qualifier.Biological_is),
    #     )
    #
    # def __eq__(self, other: Union["CVTerm", dict]) -> bool:
    #     """Compare two CVTerm objects and return boolean for equality.
    #
    #     If a dict is given, it is transformed to CVTerm.
    #     First, the qualifier is compared. If they are not identical, False is returned.
    #     Then the external resources are compared, see ExternalResources.__eq__().
    #
    #     Parameters
    #     ----------
    #     other: dict or CVTerm
    #
    #     Returns
    #     -------
    #     bool
    #         False if other is not CVTerm or dict.
    #         False if qualifiers are different, or external resources are different.
    #         True if qualifier and external resources are identical.
    #
    #     See Also
    #     --------
    #     CVTerm.from_dict()
    #     ExternalResources.__eq__()
    #     """
    #     if not isinstance(other, (CVTerm, dict)):
    #         return False
    #     if isinstance(other, dict):
    #         return self == CVTerm.from_dict(other)
    #     if self.qualifier != other.qualifier:
    #         return False
    #     if self.external_resources != other.external_resources:
    #         return False
    #     return True
    #
    # def __repr__(self) -> str:
    #     """Return the StandardizedAnnotation as str with module, class, and code to recreate it.
    #
    #     Returns
    #     -------
    #     str
    #     """
    #     return (
    #         f"{self.__class__.__module__}.{self.__class__.__qualname__}"
    #         f"({self.to_dict()})"
    #     )
    #
    # def _repr_html_(self) -> str:
    #     """Return the CVTerm as HTML string with qualifier, resources and address.
    #
    #     Returns
    #     -------
    #     str
    #         HTML formatted string
    #     """
    #     return f"""
    #                 {self.qualifier.name}:
    #                 "<p><strong>Resources</strong>"
    #                 {"<p>".join([res for res in self.external_resources.resources])}
    #                 <strong>Memory address</strong>{id(self):#x}
    #             """
    #


class StandardizedAnnotationList(UserList):
    """A list of CVTerm objects.

    Representation of multiple CVTerm objects in a list.  It is list that contains
    CVTerm objects. As a list, it means that objects can repeat, and that the order is
    maintained.

    CVTermList is built using UserList, which means that the actual list can be
    accessed using CVTermList.data. As a list, the order is kept, items may repeat.
    __init__, __setitem__, __append__, __extend__ will check the data to make sure it
    is a CVTerm or can be transformed to CVTerm using _check_CVTerm().
    All list functions that are not overloaded will behave like standard lists.

    Parameters
    ----------
    data : list
        a list containing qualifier and external resources in CVTerm format

    1. The only way to add annotation data via old format is by
       using the method "add_simple_annotations()". This will set all qualifiers as
       "bqb_is". If you want to use other qualifiers and/or nested data, use
       add_cvterms() or extend().
    2. Multiple CVTerm data can be added by using add_cvterms() or extend(). Both
       accept iterables, including CVTermList.
    """

    def __init__(
        self, data: Optional[Iterable[Union[StandardizedAnnotation, str]]] = None
    ):
        """Initialize CVTermList object.

        Parameters
        ----------
        data: Iterable of dict or CVTerm
            Dicts will be transformed to CVTerm via _check_CVTerm.

        Notes
        -----
        _check_CVTerm will raise TypeError if given a class other than dict or CVTerm,
        so initialization of CVTermList may raise TypeError.
        """
        if data is None:
            data = []
        checked_data = [
            filtered_entry
            for entry in data
            if (filtered_entry := self._check_standardized_annotation(entry))
            is not None
        ]
        super().__init__(checked_data)

    @staticmethod
    def _check_standardized_annotation(
        ann: Optional[Union[StandardizedAnnotation, str]],
    ) -> Optional["StandardizedAnnotation"]:
        if ann is None:
            return None
        if isinstance(ann, StandardizedAnnotation):
            return ann
        elif isinstance(ann, str):
            return StandardizedAnnotation(ann)
        else:
            raise TypeError(
                f"Allowed types for StandardizedAnnotationList are str and"
                f"StandardizedAnnotation, not {type(ann)}: {ann}"
            )

    @staticmethod
    def from_data(
        data: Optional[
            Union[List, Dict, "StandardizedAnnotation", "StandardizedAnnotationList"]
        ],
    ) -> "StandardizedAnnotationList":
        """Parse a CVTermList object from given data.

        Parameters
        ----------
        data: list, dict, CVTerm or CVTermList or None, optional
            This will be transformed to CVTermList class.
            None will result in an empty CVTermList.
            CVTerm and dict will be placed in a list and become a CVTermList.
            If given CVTermList, will return the data untransformed.

        Returns
        -------
        CVTermList

        Raises
        ------
        TypeError
            If not given None, dict, CVTerm or CVTermList.
        """
        if data is None:
            return StandardizedAnnotationList()
        elif isinstance(data, list):
            return StandardizedAnnotationList(data)
        elif isinstance(data, StandardizedAnnotation):
            return StandardizedAnnotationList([data])
        # elif isinstance(data, dict):
        #     return CVTermList([data])
        # TODO: Dict
        elif isinstance(data, StandardizedAnnotationList):
            return data
        else:
            raise TypeError(f"Invalid format for StandardizedAnnotationList: '{data}'")

    # def to_list_of_dicts(self) -> List[dict]:
    #     """Represent a CVTermList object as a list of python dicts.
    #
    #     Returns:
    #     -------
    #     list:
    #         a list where each item is a dict, made by CVTerm.to_dict(). Used for JSON
    #         and YAML export.
    #
    #     See Also
    #     --------
    #     CVTerm.to_dict()
    #     """
    #     return [cvterm.to_dict() for cvterm in self.data]
    #
    def add(self, ann: Iterable[Union[StandardizedAnnotation, str]]) -> None:
        """Add multiple CVTerm to CVTermList.

        Parameters
        ----------
        cvterms : Iterable
            CVTermList list of CVTerms or CVTerm dicts to be added to the CVTermList
        """
        checked_ann = [
            filtered_entry
            for entry in ann
            if (filtered_entry := self._check_standardized_annotation(entry))
            is not None
        ]
        self.extend(checked_ann)

    # def add_simple_annotations(self, data: Optional[Dict] = None) -> None:
    #     """Add simple annotation.
    #
    #     Adds standardized via old annotation format (dictionary like format).
    #     The default qualifier, i.e "bqb_is", will be used.
    #     This function will add identifiers.org to the keys given to form  the URI.
    #     If the annotation does not match the identifiers format, you should use
    #     add_cvterms directly (and create the correct link) and/or use the
    #     custompairs field of annotation.
    #
    #     This function will skip "sbo" keys since they should be added via
    #     annotation["sbo"], see MetaData.
    #
    #     Parameters
    #     ----------
    #     data : dict
    #         the data in old annotation format
    #         keys are str representing namespace
    #         If the value is a list, each value is added to the namespace, making a
    #         CVTerm with one ExternalResources object htat has mulitple URIs. If the
    #         value is a string, then only one value will be added to the namespace.
    #
    #     Examples
    #     --------
    #     >>> from cobra import Species
    #     >>> s = Species()
    #     >>> s.annotation.standardized.add_simple_annotations({"chebi": "CHEBI:17234"})
    #     >>> chebi_ent = ["CHBEI:1723456", "CHEBI:172345"]
    #     >>> s.annotation.standardized.add_simple_annotations({"chebi": chebi_ent})
    #     >>> s.annotation
    #     >>> s.annotation.standardized
    #     >>> s.annotation.annotations
    #     """
    #     if data is None:
    #         data = {}
    #
    #     if not isinstance(data, dict):
    #         raise TypeError(f"The data passed must be of type dict: {data}")
    #
    #     cvterm_list = []
    #     for key, value in data.items():
    #
    #         if not isinstance(value, (list, str)):
    #             raise TypeError(
    #                 f"The value passed must be of type list or str: {value}"
    #             )
    #         if key.lower() == "sbo":
    #             continue
    #
    #         qual = Qualifier.Biological_is
    #         # if there is only one identifier i.e. annotation
    #         # of the form { "chebi": ["CHEBI:17234"]}
    #         if isinstance(value, str):
    #             uri = ["https://identifiers.org/" + key + "/" + value]
    #         # if there are multiple identifiers for this key i.e. annotation
    #         # of the form { "chebi": ["CHEBI:124", "CHEBI:17234"]}
    #         elif isinstance(value, list):
    #             uri = [
    #                 "https://identifiers.org/" + key + "/" + identifier
    #                 for identifier in value
    #             ]
    #         else:
    #             raise TypeError(
    #                 f"The identifier passed must be of type string or list: {value}"
    #             )
    #         cvterm_list.append(CVTerm(ex_res=ExternalResources(uri), qualifier=qual))
    #     self.add_cvterms(cvterm_list)
    #
    # def delete_annotation(self, resource: Union[str, Pattern]) -> None:
    #     r"""Delete annotation - the converse of add_simple_annotation.
    #
    #     This will go over the CVTerms, and delete all resources that match the pattern.
    #     It will call the funciton recursively for ExternalResources that have
    #     nested_data.
    #     CVTerms that end with neither resources nor nested data are removed.
    #
    #     Parameters
    #     ----------
    #     resource: str or Pattern
    #
    #     Examples
    #     --------
    #     >>> from cobra.io import load_model
    #     >>> e_coli = load_model('iJO1366')
    #     >>> e_coli.annotation
    #     >>> e_coli.annotation.standardized.delete_annotation('bigg')
    #     >>> e_coli.annotation
    #     >>> e_coli.annotation.annotations
    #     >>> e_coli.metabolites[0].annotation.standardized
    #     >>> e_coli.metabolites[0].annotation.standardized.delete_annotation(r'/bi\S+')
    #     """
    #     regex_searcher = re.compile(resource)
    #     tmp_cvterm_list = []
    #     for cvterm in self.data:
    #         cvterm.external_resources.resources = [
    #             res
    #             for res in cvterm.external_resources.resources
    #             if not regex_searcher.findall(res)
    #         ]
    #         if cvterm.external_resources.nested_data:
    #             cvterm.external_resources.nested_data.delete_annotation(resource)
    #         if (
    #             cvterm.external_resources.resources
    #             or cvterm.external_resources.nested_data
    #         ):
    #             tmp_cvterm_list.append(cvterm)
    #     self.data = tmp_cvterm_list
    #
    # @property
    # def annotations(self) -> Dict:
    #     """Return CVTermList as annotation dictionary.
    #
    #     This function will return the CVTermList as a sorted annotation dictionary in
    #     the older annotation format. In the dictionary, the keys are the namespaces,
    #     while the values are lists of annotation identifiers.
    #
    #     Qualifiers are not present in the annotation dictionary. This function will
    #     use CVTermList.resources(), which will get all resources of all CVTerm objects,
    #     including nested resources.
    #
    #     Different CVTerms will be unified by namespace. Namespaces and values are
    #     sorted in the dictionary, to avoid shifting caused by usage of sets.
    #
    #     For example, a CVTermList that looks like
    #
    #     [
    #         {
    #             "external_resources": {
    #                 "resources": [
    #                     "https://identifiers.org/uniprot/P69906",
    #                     "https://identifiers.org/uniprot/P68871",
    #                     "https://identifiers.org/kegg.compound/C00032",
    #                 ]
    #             },
    #             "qualifier": "bqb_hasPart",
    #         },
    #         {
    #             "qualifier": "bqb_hasPart",
    #             "external_resources": {
    #                 "resources": [
    #                     "https://identifiers.org/uniprot/P69905",
    #                     "https://www.uniprot.org/uniprot/P68871",
    #                     "https://identifiers.org/chebi/CHEBI:17627",
    #                 ],
    #             "nested_data": {
    #                 "qualifier": "bqb_isDescribedBy",
    #                 "external_resources": {
    #                     "resources": [
    #                         "https://identifiers.org/eco/000000",
    #                     ]
    #                 },
    #             },
    #         },
    #     ]
    #
    #     Will be outputted as a dictionary that looks like
    #     {
    #         "chebi": ["CHEBI:17627"],
    #         "eco": ["000000"],
    #         "kegg.compound": ["C00032"],
    #         "uniprot": ["P68871", "P69905", "P69906"],
    #     }
    #
    #     Returns
    #     -------
    #     dict
    #         Dictionary where keys are namespaces, sorted in ascending order. Values
    #         are lists of identifiers, also sorted in ascentding order.
    #
    #     """
    #     annotation_dict = {}
    #     resources = self.resources
    #     for res in resources:
    #         if re.match(URL_IDENTIFIERS_PATTERN, res):
    #             identifier_match = parse_identifiers_uri(res)
    #             if identifier_match is None:
    #                 continue
    #             namespace, identifier = identifier_match
    #             if namespace in annotation_dict.keys():
    #                 annotation_dict[namespace].append(identifier)
    #             else:
    #                 annotation_dict[namespace] = [identifier]
    #     return {k: sorted(annotation_dict[k]) for k in sorted(annotation_dict.keys())}
    #
    # @property
    # def resources(self) -> FrozenSet[str]:
    #     """Get all resources.
    #
    #     Returns:
    #     -------
    #     FrozenSet:
    #         a set of all external resources in the original self.data list of CVTerms
    #         including external resources of nested data. The Set contains the URIs as
    #         strings, not in the ExternalResources format.
    #     """
    #     resources = set()
    #     for datum in self.data:
    #         if (
    #             datum.external_resources.resources
    #             or datum.external_resources.nested_data
    #         ):
    #             resources.update(datum.external_resources.resource_nested)
    #     return frozenset(resources)
    #
    # @property
    # def qualifiers(self) -> FrozenSet[Qualifier]:
    #     """Get all qualifiers used by CVTerm objects in the CVTermList.
    #
    #     Returns:
    #     -------
    #     FrozenSet:
    #         a frozen set of all qualifiers in the original self.data list of CVTerms
    #     """
    #     qualifier_set = set()
    #     for datum in self.data:
    #         qualifier_set.add(datum.qualifier)
    #     return frozenset(qualifier_set)
    #
    # def query(
    #     self,
    #     search_function: Union[str, Pattern, Callable],
    #     attribute: Union[str, None] = None,
    # ) -> "CVTermList":
    #     """Query the CVTermList and return a list of CVTerm objects.
    #
    #     Parameters
    #     ----------
    #     search_function : a string, regular expression or function
    #         Used to find the matching elements in the list.
    #         - a regular expression (possibly compiled), in which case the
    #         given attribute of the object should match the regular expression.
    #         - a function which takes one argument and returns True for
    #         desired values
    #
    #     attribute : string or None
    #         the name attribute of the object to passed as argument to the
    #         `search_function`. If this is None and a regular expression/string is given,
    #          will match the regular expression to both qualifier and resources.
    #
    #     Returns
    #     -------
    #     CVTermList
    #         a new list of CVTerm objects which match the query
    #
    #     Examples
    #     --------
    #     >>> from cobra.io import load_model
    #     >>> model = load_model('iJO1366')
    #     >>> model.annotation.standardized.query('bqb', 'qualifier')
    #     >>> import re
    #     >>> regex = re.compile('^bqm')
    #     >>> model.annotation.standardized.query(regex, 'qualifier')
    #     """
    #
    #     def select_attribute(
    #         x: CVTerm,
    #     ) -> Union[CVTerm, ExternalResources, Qualifier, set]:
    #         if attribute is None:
    #             return x
    #         else:
    #             return getattr(x, attribute)
    #
    #     try:
    #         # if the search_function is a regular expression
    #         regex_searcher = re.compile(search_function)
    #         print(f"Search function: '{search_function}'")
    #         if attribute is None:
    #             attribute = ""
    #
    #         if attribute == "qualifier":
    #             matches = [
    #                 cvterm
    #                 for cvterm in self.data
    #                 if (
    #                     regex_searcher.findall(select_attribute(cvterm).name) != []
    #                     or regex_searcher.findall(select_attribute(cvterm).value) != []
    #                 )
    #             ]
    #         elif attribute == "resources":
    #             matches = [
    #                 cvterm
    #                 for cvterm in self.data
    #                 if any(
    #                     regex_searcher.findall(res) for res in select_attribute(cvterm)
    #                 )
    #             ]
    #         elif attribute == "external_resources":
    #             matches = [
    #                 cvterm
    #                 for cvterm in self.data
    #                 if any(
    #                     regex_searcher.findall(res)
    #                     for res in select_attribute(cvterm).resources
    #                 )
    #             ]
    #         else:
    #             matches = [
    #                 cvterm
    #                 for cvterm in self.data
    #                 if regex_searcher.findall(cvterm.qualifier.name) != []
    #                 or regex_searcher.findall(cvterm.qualifier.value) != []
    #                 or any(regex_searcher.findall(res) for res in cvterm.resources)
    #             ]
    #     except TypeError:
    #         matches = [
    #             cvterm
    #             for cvterm in self.data
    #             if search_function(select_attribute(cvterm))
    #         ]
    #
    #     results = self.__class__(matches)
    #     return results
    #
    def __setitem__(self, key: int, value: Union[StandardizedAnnotation, str]) -> None:
        """Set item in CVTermList.

        This function will check that it is a valid item, converting it to CVTerm if
        necessary (see _check_CVTerm).

        Parameters
        ----------
        key: int
        value: CVTerm or dict
            dict will be converted to CVTerm
        """
        checked_value = self._check_standardized_annotation(value)
        if checked_value is None:
            TypeError(f"Value cannot be None.")
            # TODO: Elaborate (or automatically delete when None)
        UserList.__setitem__(self, key, checked_value)

    def append(self, item: Union[StandardizedAnnotation, str]) -> None:
        """Append CVTerm to end.

        Parameters
        ----------
        item: CVTerm or dict
            dict will be converted to CVTerm.
        """
        checked_item = self._check_standardized_annotation(item)
        if checked_item is None:
            TypeError(f"Item cannot be None.")
            # TODO: Elaborate (or do nothing when None)
        UserList.append(self, checked_item)

    def extend(
        self,
        iterable: Union[
            "StandardizedAnnotationList", Iterable[Union[StandardizedAnnotation, str]]
        ],
    ) -> None:
        """Extend data list by appending elements from the iterable.

        Parameters
        ----------
        iterable : Iterable
        """
        if isinstance(iterable, StandardizedAnnotationList):
            self.data.extend(iterable.data)
        elif isinstance(iterable, Iterable):
            self.data.extend(
                [
                    checked_item
                    for item in iterable
                    if (checked_item := self._check_standardized_annotation(item))
                    is not None
                ]
            )

    # def __eq__(self, other: Union[list, "CVTermList"]) -> bool:
    #     """Compare two CVTermList objects to find out whether they are the same.
    #
    #     Equality is defined as them having the same data, but not necessarily the same
    #     objects. If the given item is not a CVTermList or list, this function will
    #     return False.
    #
    #     Parameters
    #     ----------
    #     other: CVTermList or list
    #
    #     Returns
    #     -------
    #     bool: True if the data matches, False otherwise
    #     """
    #     if isinstance(other, list):
    #         return self.__eq__(CVTermList.from_data(other))
    #     if not isinstance(other, CVTermList):
    #         return False
    #     if len(self.data) != len(other.data):
    #         return False
    #     for self_i, other_i in zip(self.data, other.data):
    #         if self_i != other_i:
    #             return False
    #     return True
    #
    # def _repr_html_(self) -> str:
    #     """Generate CVTermList as HTML.
    #
    #     Returns
    #     -------
    #     str
    #         HTML representation of the list of CVTerm resources.
    #     """
    #     return f"""CVTermList{"<p>".join([cvterm._repr_html_()
    #                                       for cvterm in self.data])}"""
    #


# class ExternalResources:
#     """Representation of a single set of resources and its nested annotation.
#
#     Each resource in the resources fields is a URI that uniquely identifies both the
#     resource and the data within the resource. Since a URI is not a URL, it does not
#     have to map to a physical Web object; it simply needs to identify, uniquely, a
#     controlled vocabulary term or database object.
#     These URIs are MIRIAM identifiers, following the format defined in
#     https://identifiers.org/.
#
#     The format allowed is
#     https://identifiers.org/[provider_code/]namespace:accession
#
#     The optional parameter provider_code denotes the Provider Code part of the Prefix.
#     It is trailed by a slash, to separate it from the required namespace,
#      which is followed by a colon and the accession.
#
#     Example Resources URIs
#     ----------------------
#     https://identifiers.org/pubmed:22140103
#     https://identifiers.org/ec-code:1.1.1.1
#
#     https://identifiers.org/epmc/pubmed:22140103
#     https://identifiers.org/expasy/ec-code:1.1.1.1
#
#
#     Parameters
#     ----------
#     resources: list or str
#         A list of URIs (resources)
#     nested_data : CVTermList
#         Nested annotation, in CVTermList format
#
#     Can also be created from dictionary, see ExternalResources.from_dict()
#     """
#
#     def __init__(
#         self,
#         resources: Optional[Union[List, str]] = None,
#         nested_data: Optional[Union[Dict, CVTerm, CVTermList]] = None,
#     ):
#         """Initialize ExternalResources object.
#
#         Parameters
#         ----------
#         resources: list or str, optional
#             A list of URIs (resources) or str. Str will be placed in a list if only
#             one resource is used.
#         nested_data : dict or CVTerm or CVTermList, optional
#             Nested annotation, in CVTermList format
#
#         See Also
#         --------
#         ExternalResources.resources()
#         ExternalResources.nested_data()
#         """
#         self._resources = None
#         self._nested_data = None
#         self.resources = resources
#         if resources:
#             self.nested_data = nested_data
#
#     @property
#     def resources(self) -> FrozenSet:
#         """Get resources of ExternalResources.
#
#         Returns
#         -------
#         frozenset:
#             The list of URIs in a frozenset.
#         """
#         if self._resources is None:
#             return frozenset()
#         return frozenset(self._resources)
#
#     @resources.setter
#     def resources(self, value: Optional[Union[List[str], str]]) -> None:
#         """Set resources of ExternalResources.
#
#         Will set the URIs of ExternalResources.
#
#         Parameters
#         ----------
#         value: str or list
#             Should be list of URIs. If only one string is given, the function
#             assumes it is one URI, and it is placed in a list.
#             If None is given, the resources field is set to an empty list.
#
#         Raises
#         ------
#         TypeError
#             If value is neither list nor str.
#         """
#         if value is None:
#             value = []
#         if isinstance(value, list):
#             self._resources = value
#         elif isinstance(value, str):
#             self._resources = [value]
#         else:
#             raise TypeError(f"The resources must be a string or a list: {value}")
#
#     @property
#     def resource_nested(self) -> FrozenSet:
#         """Get resources, including resources of nested data.
#
#         This will get resources of the current ExternalResource and all resources
#         of nested data, if they exist.
#
#         Returns
#         -------
#         FrozenSet: a frozen set of URIs
#         """
#         resources = set()
#         resources.update(self.resources)
#         if self.nested_data:
#             resources.update(self.nested_data.resources)
#         return frozenset(resources)
#
#     @property
#     def nested_data(self) -> Optional[CVTermList]:
#         """Get nested data of ExternalResources.
#
#         Returns
#         -------
#         CVTermList: optional
#             Will return None if this ExternalResources object has no nested data.
#         """
#         return self._nested_data
#
#     @nested_data.setter
#     def nested_data(self, value: Optional[Union[Dict, List, CVTerm, CVTermList]]):
#         """Set nested data of ExternalResources.
#
#         This function will convert the given value to CVTermList as follows
#         - dict is converted to CVTerm via CVTerm.from_dict()
#         - list is converted via CVTermList.from_data()
#         - CVTerm is placed in a list and tranformed to CVTermList
#         - CVTermList is unchanged
#         - None will set the nested_data as an empty CVTermList (used for new
#             ExternalResources objects)
#
#         Parameters
#         ----------
#         value: dict or list or CVTerm or CVTermList, optional
#
#         Raises
#         ------
#         TypeError
#             If value is neither None, list, CVTerm or CVTermList.
#         """
#         if value is None:
#             self._nested_data = CVTermList()
#         elif isinstance(value, CVTerm):
#             self._nested_data = CVTermList([value])
#         elif isinstance(value, CVTermList):
#             self._nested_data = value
#         elif isinstance(value, list):
#             self._nested_data = CVTermList.from_data(value)
#         elif isinstance(value, dict):
#             self._nested_data = CVTermList.from_data(CVTerm.from_dict(value))
#         else:
#             raise TypeError(
#                 f"The nested data structure does not have valid CVTerm format: {value}"
#             )
#
#     def to_dict(self) -> Dict[str, Any]:
#         """Generate a dict representing an ExternalResource object.
#
#         Returns
#         -------
#         dict
#             A dictionary with two keys
#             "resources" - sorted list of resources.
#             "nested_data" - optional, dict of nested_data. See CVTermList.to_dict()
#         """
#         if self._resources is None:
#             return {}
#         ex_dic: Dict[str, Any] = {"resources": sorted(self._resources)}
#         if self.nested_data is not None and len(self.nested_data):
#             ex_dic["nested_data"] = self.nested_data.to_list_of_dicts()
#         return ex_dic
#
#     @classmethod
#     def from_dict(cls, input_data: Dict) -> "ExternalResources":
#         """Generate ExternalResources object from dict.
#
#         Parameters
#         ----------
#         input_data: dict
#             Accepts a dictionary as input that has two keys
#             "resources", optional a list of URIs
#             "nested_data", optional - nested data. Nested data must be a list,
#             CVTerm or CVTermList.
#
#         Returns
#         -------
#         ExternalResources
#
#         See Also
#         --------
#         ExternalResources.nested_data
#
#         """
#         ex_res = cls(
#             resources=input_data.get("resources", None),
#             nested_data=input_data.get("nested_data", None),
#         )
#         return ex_res
#
#     def __eq__(self, other: "ExternalResources") -> bool:
#         """Compare two ExternalResources object.
#
#         Compare two ExternalResources objects to find out whether
#         they are same (have same data) or not.
#
#         Returns
#         -------
#         bool
#             True if resources are equal in both objects, and both objects have the same
#             nested_data (None or identical CVTerm). False otherwise.
#         """
#         if self.resources != other.resources:
#             return False
#         if self.nested_data is None and other.nested_data is None:
#             return True
#         elif self.nested_data is None or other.nested_data is None:
#             return False
#         elif self.nested_data != other.nested_data:
#             return False
#         return True
#
#     def __str__(self) -> str:
#         """Return ExternalResources dictionary as str.
#
#         Returns
#         -------
#         str
#             String output of to_dict().
#
#         See Also
#         --------
#         ExternalResources.to_dict()
#         """
#         return str(self.to_dict())
#
#     def __repr__(self) -> str:
#         """Generate ExternalResources as str, with class.
#
#         Returns
#         -------
#         str
#         """
#         return (
#             f"{self.__class__.__module__}.{self.__class__.__qualname__} "
#             f"{self.__str__()}"
#         )
#
#     def _repr_html_(self) -> str:
#         """Generate ExternalResources as HTML.
#
#         Returns
#         -------
#         str
#             HTML representation of the external resource.
#         """
#         html_rep = f"""
#                     <p>
#                         <strong>Resources</strong>{"<p>".join(self.resources)}<p>"""
#         if self.nested_data is not None and len(self.nested_data):
#             html_rep += f"""<strong>Nested Data</strong>{self.nested_data}<p>"""
#
#         html_rep += f"""
#                         <strong>Memory address</strong>{id(self):#x}
#                 <p>"""
#         return html_rep
