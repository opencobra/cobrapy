"""Define the Controlled Vocabulary term class."""

import re
from collections import UserList
from collections.abc import Iterable as ABCIterable, KeysView, MutableMapping
from typing import (
    Any,
    Callable,
    Dict,
    FrozenSet,
    Iterable,
    List,
    Optional,
    Pattern,
    Tuple,
    Union,
)

from cobra.core.metadata.identifier import get_default_qualifier

from .. import object as CObject
from cobra.core.metadata import Identifier, Qualifier, Identifier


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
        self._parent = None
        # TODO: Keep track of target Object, so we can do annotation.remove

    def _set_parent(self, parent):
        if self._parent is None or self._parent is parent:
            self._parent = parent
        else:
            raise ValueError(
                "StandardizedAnnotation already has a different parent. Create a new "
                "annotation if you would like to add an annotation to a second object."
            )

    def remove_from_parent(self):
        if self._parent is None:
            raise ValueError(
                "Cannot remove annotation, since no object is associated with annotation."
            )
        self._parent.remove(self)
        self._parent = None

    def _remove_identifier(self, identifier):
        self.identifiers.remove(identifier)

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
        for idf in self._identifiers:
            idf._set_parent(self)

    def add_identifiers(self, identifiers: Iterable[Union[str, "Identifier"]]) -> None:
        identifiers = self.check_identifier_type(identifiers)
        for idf in identifiers:
            idf._set_parent(self)
        self.identifiers.extend(identifiers)

    @property
    def uris(self) -> FrozenSet[str]:
        l = {entry.uri for entry in self.identifiers}
        for entry in self.annotations:
            l.update(entry.uris)
        return frozenset(l)

    @property
    def annotations(self) -> "StandardizedAnnotationList":
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
            Union[
                "Identifier",
                str,
                Dict[str, str],
                Tuple[str, str],
                Iterable[Union["Identifier", str, Dict[str, str], Tuple[str, str]]],
            ]
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
        elif isinstance(identifiers, (Identifier, str, dict)):
            return [Identifier.from_data(identifiers)]
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
    ) -> "StandardizedAnnotationList":
        return StandardizedAnnotationList.from_data(ann)
        #
        # if ann is None:
        #     return StandardizedAnnotationList.from_data(None)
        # elif isinstance(ann, StandardizedAnnotation):
        #     return [ann]
        # elif isinstance(ann, ABCIterable):
        #     return [x for y in ann for x in __class__.check_annotation_type(y)]
        # else:
        # raise TypeError(
        #     f"Allowed types for nested StandardizedAnnotation annotations"
        #     f"are StandardizedAnnotation or a list of StandardizedAnnotation"
        #     f"objects, not {type(ann)}: ann"
        # )
        #

    def to_dict(self) -> Dict:
        """Represent a CVTerm object in python dict.

        Returns
        -------
        dict:
            A dict that has two keys
            "qualifier" - the qualifier as a string
            "external_resources" - the resources as a dictionary

        See Also
        --------
        ExternalResources.to_dict()

        """
        d = {
            "qualifier": self.qualifier.value,
            "identifiers": [identifier.to_dict() for identifier in self.identifiers],
        }
        if self.annotations:
            d["annotations"] = self.annotations.to_list_of_dicts()
        return d

    def to_tuples(self) -> List[Tuple[str, str]]:
        return [(idf.namespace, idf.identifier) for idf in self.identifiers]

    def to_records(self) -> List[Dict[str, Any]]:
        l, _ = self._to_records()
        return l

    def _to_records(
        self, group_counter: int = 1, parent_group: int = 0
    ) -> Tuple[List[Dict[str, Any]], int]:
        l = []
        for identifier in self.identifiers:
            entry = {
                "qualifier": self.qualifier.value,
                "uri": identifier.uri,
                "namespace": identifier.namespace,
                "identifier": identifier.identifier,
                "annotation_group": group_counter,
                "parent_group": parent_group,
            }
            l.append(entry)
        if self.annotations:
            new_l, group_counter = self.annotations._to_records(
                group_counter=(group_counter + 1), parent_group=group_counter
            )
            l.extend(new_l)
        else:
            group_counter = group_counter + 1
        return l, group_counter

    @classmethod
    def from_dict(cls, data_dict: Dict) -> "StandardizedAnnotation":
        """Generate a CVTerm object based on a python dict.

        Parameters
        ----------
        data_dict: dict
            A dict that has two keys
            "qualifier" - the qualifier as a string, optional. If not present, the
            qualifier is set to bqb_is.
            "external_resources" - the resources as a dictionary, optional

        Returns
        -------
        CVTerm

        See Also
        --------
        ExternalResources.to_dict()

        """
        return cls(
            identifiers=data_dict.get("identifiers", None),
            qualifier=data_dict.get("qualifier", Qualifier.Biological_is),
            annotations=data_dict.get("annotations", None),
        )

    def __eq__(self, other: Any) -> bool:
        """Compare two CVTerm objects and return boolean for equality.

        If a dict is given, it is transformed to CVTerm.
        First, the qualifier is compared. If they are not identical, False is returned.
        Then the external resources are compared, see ExternalResources.__eq__().

        Parameters
        ----------
        other: dict or CVTerm

        Returns
        -------
        bool
            False if other is not CVTerm or dict.
            False if qualifiers are different, or external resources are different.
            True if qualifier and external resources are identical.

        See Also
        --------
        CVTerm.from_dict()
        """

        if isinstance(other, dict):
            return self == StandardizedAnnotation.from_dict(other)
        if isinstance(other, StandardizedAnnotation):
            if self.qualifier != other.qualifier:
                return False
            if len(self.identifiers) != len(other.identifiers):
                return False
            for idf in self.identifiers:
                if not idf in other.identifiers:
                    return False
            if self.annotations != other.annotations:
                return False
            return True
        return False

    def __repr__(self) -> str:
        """Return the StandardizedAnnotation as str with module, class, and code to recreate it.

        Returns
        -------
        str
        """
        return (
            f"{self.__class__.__module__}.{self.__class__.__qualname__}"
            f"({self.to_dict()})"
        )

    def _repr_html_(self) -> str:
        """Return the CVTerm as HTML string with qualifier, resources and address.

        Returns
        -------
        str
            HTML formatted string
        """
        # TODO: Fix this HTML
        return f"""
                    {self.qualifier.name}:
                    <p><strong>Identifiers</strong>
                    {"</p><p>".join([res._repr_html() for res in self.identifiers])}
                    <p><strong>Annotations</strong></p>
                    <p>{self.annotations._repr_html_()}</p>
                    <strong>Memory address</strong>{id(self):#x}
                """


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
        self,
        data: Optional[
            Iterable[Union[StandardizedAnnotation, Dict, Identifier, str]]
        ] = None,
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
            if (not isinstance(entry, (str, Identifier)))
            and (filtered_entry := self._check_standardized_annotation(entry))
            is not None
        ]
        # str and Identifier instances are handled separately and added as one
        # Standardized annotation instance with the default qualifier (Biological_is).
        if no_qualifier_data := [
            entry for entry in data if isinstance(entry, (str, Identifier))
        ]:
            checked_data.insert(
                0, StandardizedAnnotation(identifiers=no_qualifier_data)
            )
        for entry in checked_data:
            entry._set_parent(self)
        super().__init__(checked_data)

    @staticmethod
    def _check_standardized_annotation(
        ann: Optional[Union[StandardizedAnnotation, Dict, str]],
    ) -> Optional["StandardizedAnnotation"]:
        if ann is None:
            return None
        if isinstance(ann, StandardizedAnnotation):
            return ann
        elif isinstance(ann, str):
            return StandardizedAnnotation(ann)
        elif isinstance(ann, dict):
            return StandardizedAnnotation.from_dict(ann)
        else:
            raise TypeError(
                f"Allowed types for StandardizedAnnotation are str and"
                f"StandardizedAnnotation, not {type(ann)}: {ann}"
            )
        # TODO: Handle dict

    @staticmethod
    def from_data(
        data: Optional[
            Union[
                Iterable[Union[str, Dict, "StandardizedAnnotation"]],
                str,
                Dict,
                "StandardizedAnnotation",
                "StandardizedAnnotationList",
            ]
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
        elif isinstance(data, StandardizedAnnotationList):
            return data
        elif isinstance(data, (StandardizedAnnotation, dict, str)):
            return StandardizedAnnotationList([data])
        elif isinstance(data, ABCIterable):
            return StandardizedAnnotationList(data)
        else:
            raise TypeError(f"Invalid format for StandardizedAnnotationList: '{data}'")

    def to_list_of_dicts(self) -> List[dict]:
        """Represent a CVTermList object as a list of python dicts.

        Returns:
        -------
        list:
            a list where each item is a dict, made by CVTerm.to_dict(). Used for JSON
            and YAML export.

        See Also
        --------
        CVTerm.to_dict()
        """
        return [cvterm.to_dict() for cvterm in self.data]

    def to_records(self):
        l, _ = self._to_records()
        return l

    def _to_records(
        self, group_counter: int = 1, parent_group: int = 0
    ) -> Tuple[List[Dict], int]:
        l = []
        for entry in self.data:
            new_l, group_counter = entry._to_records(
                group_counter=group_counter, parent_group=parent_group
            )
            l.extend(new_l)
        return l, group_counter

    def _find_first_by_qualifier(
        self,
        qualifier: Qualifier = Qualifier.Biological_is,
    ) -> Optional[StandardizedAnnotation]:
        for entry in self.data:
            if entry.qualifier == qualifier:
                return entry
        return None

    def _find_first_or_create_by_qualifier(
        self, qualifier: Qualifier = Qualifier.Biological_is
    ) -> StandardizedAnnotation:
        entry = self._find_first_by_qualifier(qualifier)
        if entry is None:
            entry = StandardizedAnnotation(qualifier=qualifier, identifiers=[])
            entry._set_parent(self)
            self.data.insert(0, entry)
        return entry

    def add(self, ann: Iterable[Union[StandardizedAnnotation, Dict, str]]) -> None:
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
    @property
    def identifiers(self) -> FrozenSet[Identifier]:
        """Get all Identifiers.

        Returns:
        -------
        FrozenSet:
            a set of all external resources in the original self.data list of CVTerms
            including external resources of nested data. The Set contains the URIs as
            strings, not in the ExternalResources format.
        """
        resources = set()
        for entry in self.data:
            resources.update(entry.identifiers)
            if entry.annotations:
                resources.update(entry.annotations.identifiers)
        return frozenset(resources)

    @property
    def uris(self) -> FrozenSet[str]:
        l = set()
        for entry in self.data:
            l.update(entry.uris)
        return frozenset(l)

    @property
    def qualifiers(self) -> FrozenSet[Qualifier]:
        """Get all qualifiers used by CVTerm objects in the CVTermList.

        Note it does not return nested qualifiers.

        Returns:
        -------
        FrozenSet:
            a frozen set of all qualifiers in the original self.data list of CVTerms
        """
        qualifier_set = set()
        for entry in self.data:
            qualifier_set.add(entry.qualifier)
        return frozenset(qualifier_set)

    def __iter__(self):
        return iter(self.data)

    def __len__(self):
        return len(self.data)

    def query(
        self,
        search_function: Union[str, Pattern, Callable],
        attribute: Union[str, None] = None,
    ) -> "StandardizedAnnotationList":
        """Query the CVTermList and return a list of CVTerm objects.

        Parameters
        ----------
        search_function : a string, regular expression or function
            Used to find the matching elements in the list.
            - a regular expression (possibly compiled), in which case the
            given attribute of the object should match the regular expression.
            - a function which takes one argument and returns True for
            desired values

        attribute : string or None
            the name attribute of the object to passed as argument to the
            `search_function`. If this is None and a regular expression/string is given,
             will match the regular expression to both qualifier and resources.

        Returns
        -------
        CVTermList
            a new list of CVTerm objects which match the query

        Examples
        --------
        >>> from cobra.io import load_model
        >>> model = load_model('iJO1366')
        >>> model.annotation.standardized.query('bqb', 'qualifier')
        >>> import re
        >>> regex = re.compile('^bqm')
        >>> model.annotation.standardized.query(regex, 'qualifier')
        """

        # TODO: Clean up this whole method.
        def select_attribute(
            x: StandardizedAnnotation,
        ) -> Union[StandardizedAnnotation, Identifier, Qualifier, set]:
            if attribute is None:
                return x
            else:
                return getattr(x, attribute)

        try:
            # if the search_function is a regular expression
            regex_searcher = re.compile(search_function)
            print(f"Search function: '{search_function}'")
            if attribute is None:
                attribute = ""

            if attribute == "qualifier":
                matches = [
                    cvterm
                    for cvterm in self.data
                    if (
                        regex_searcher.findall(select_attribute(cvterm).name) != []
                        or regex_searcher.findall(select_attribute(cvterm).value) != []
                    )
                ]
            elif attribute == "identifiers":
                matches = [
                    cvterm
                    for cvterm in self.data
                    if any(
                        regex_searcher.findall(res.uri)
                        for res in select_attribute(cvterm)
                    )
                ]
            else:
                matches = [
                    cvterm
                    for cvterm in self.data
                    if regex_searcher.findall(cvterm.qualifier.name) != []
                    or regex_searcher.findall(cvterm.qualifier.value) != []
                    or any(
                        regex_searcher.findall(res.uri) for res in cvterm.identifiers
                    )
                ]
        except TypeError as err:
            print(err)
            matches = [
                cvterm
                for cvterm in self.data
                if search_function(select_attribute(cvterm))
            ]

        results = self.__class__(matches)
        return results

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
            raise TypeError(f"Value cannot be None.")
            # TODO: Elaborate (or automatically delete when None)
        checked_value._set_parent(self)
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
            raise TypeError(f"Item cannot be None.")
            # TODO: Elaborate (or do nothing when None)
        checked_item._set_parent(self)
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
            self.extend(iterable.data)
        elif isinstance(iterable, Iterable):
            checked_data = [
                checked_item
                for item in iterable
                if (checked_item := self._check_standardized_annotation(item))
                is not None
            ]
            for d in checked_data:
                d._set_parent(self)
            self.data.extend(checked_data)

    def __eq__(self, other: Union[Iterable, "StandardizedAnnotationList"]) -> bool:
        """Compare two CVTermList objects to find out whether they are the same.

        Equality is defined as them having the same data, but not necessarily the same
        objects. If the given item is not a CVTermList or list, this function will
        return False.

        Parameters
        ----------
        other: CVTermList or list

        Returns
        -------
        bool: True if the data matches, False otherwise
        """
        if isinstance(other, ABCIterable) and not isinstance(
            other, StandardizedAnnotationList
        ):
            return self.__eq__(StandardizedAnnotationList.from_data(other))
        if not isinstance(other, StandardizedAnnotationList):
            return False
        if len(self.data) != len(other.data):
            return False
        for other_entry in other.data:
            if not other_entry in self.data:
                return False
        return True

    def _repr_html_(self) -> str:
        """Generate CVTermList as HTML.

        Returns
        -------
        str
            HTML representation of the list of CVTerm resources.
        """
        entries = [cvterm._repr_html_() for cvterm in self.data]
        return f"""StandardizedAnnotationList{"<p>".join(entries)}"""


class SimplifiedAnnotationInterface(MutableMapping):
    def __init__(self, standardized_annotations: StandardizedAnnotationList):
        self._annotations = standardized_annotations

    def add(
        self,
        data: Optional[
            Union[
                Dict,
                str,
                Tuple[str, str],
                List[Union[str, Identifier, Tuple[str, str]]],
            ]
        ] = None,
    ) -> None:
        if data is None:
            data = []

        if isinstance(data, dict):
            data = list(data.items())

        if isinstance(data, (str, tuple)):
            data = [data]

        if not isinstance(data, list):
            raise TypeError(
                "The supplied annotations were not of type List, or could "
                "not be converted to a list."
            )

        cvterms = {}
        for entry in data:
            if isinstance(entry, tuple):
                if len(entry) != 2:
                    raise ValueError(
                        f"Only tuples of length 2 can be converted to annotations."
                    )
                l = []
                if isinstance(entry[1], list):
                    l.extend([(entry[0], v) for v in entry[1]])
                else:
                    l.append(entry)
                entry = l
            elif isinstance(entry, str):
                entry = [entry]
            else:
                raise TypeError("Entry could could not be converted to an Identifier.")
            for x in entry:
                x = Identifier.from_data(x)
                if x.namespace is None:
                    raise ValueError(
                        f"Could not determine namespace of identifier {x}."
                    )

                qualifier = get_default_qualifier(x.namespace)
                if not qualifier in cvterms:
                    cvterms[qualifier] = []
                cvterms[qualifier].append(x)

        for qualifier, identifiers in cvterms.items():
            ann = self._annotations._find_first_or_create_by_qualifier(qualifier)
            ann.add_identifiers(identifiers)

    def __getitem__(
        self, idx: Union[int, str, Qualifier]
    ) -> Optional[Union[StandardizedAnnotation, Dict[str, str], List[str]]]:
        if isinstance(idx, int):
            return self._annotations.data[idx]
        if isinstance(idx, Qualifier):
            idfs = self._annotations._find_first_by_qualifier(idx)
            if idfs is None:
                return None
            idfs = idfs.to_tuples()
            d = {}
            for k, v in idfs:
                if not k in d:
                    d[k] = [v]
                else:
                    d[k].append(v)
            return d
        if isinstance(idx, str):
            qual = get_default_qualifier(idx)
            ann = self._annotations._find_first_by_qualifier(qual)
            if ann is None:
                return ann
            return [
                v
                for idf in ann.identifiers
                if (v := idf.identifier) is not None and idf.namespace == idx
            ]
        raise TypeError("Index should be of type int, str or Qualifier.")

    def __delitem__(self, idx: str):
        if not isinstance(idx, str):
            raise TypeError("Index should be of type str.")
        qual = get_default_qualifier(idx)
        ann = self._annotations._find_first_by_qualifier(qual)
        if ann is None:
            raise IndexError(
                f"Could not find annotations for f'{idx}' (qualifier '{qual}')"
            )

        for idf in list(ann.identifiers):
            if idf.namespace is None:
                continue
            if idf.namespace == idx:
                idf.remove_from_parent()

    def __setitem__(
        self, key: str, value
    ) -> Optional[Union[StandardizedAnnotation, Dict[str, str], List[str]]]:
        if not isinstance(key, str):
            raise TypeError("Index should be of type str.")

        try:
            del self[key]
        except IndexError:
            pass

        self.add({key: value})

    def __eq__(self, other: Union[Iterable, "StandardizedAnnotationList"]) -> bool:
        return self._annotations == other

    def items(self):
        visited_qualifiers = set()
        visited_namespaces = set()
        for entry in self._annotations.data:
            qualifier = entry.qualifier
            # Only the first occurence of each qualifier is handled in simplified
            # annotations.
            if qualifier in visited_qualifiers:
                continue
            visited_qualifiers.add(qualifier)

            identifiers = True
            while identifiers:
                identifiers = []
                current_namespace = None
                for identifier in entry.identifiers:
                    if identifier.namespace is None:
                        continue
                    if identifier.namespace in visited_namespaces:
                        continue
                    if current_namespace is None:
                        if qualifier == get_default_qualifier(identifier.namespace):
                            current_namespace = identifier.namespace
                            identifiers.append(identifier)
                    else:
                        if identifier.namespace == current_namespace:
                            identifiers.append(identifier)
                visited_namespaces.add(current_namespace)
                if identifiers:
                    yield (current_namespace, identifiers)
                else:
                    break

    def __iter__(self):
        for k, _ in self.items():
            yield k

    def keys(self):
        return KeysView(self)
        # for k in self:
        #     yield k

    def values(self):
        for _, v in self.items():
            yield v

    def identifiers(self):
        for _, v in self.items():
            for identifier in v:
                yield identifier

    def tuples(self):
        for k, v in self.items():
            for identifier in v:
                yield (k, identifier)

    @property
    def number_of_identifiers(self):
        return sum(1 for _ in self.identifiers())

    def __len__(self):
        return sum(1 for _ in self)

    def delete_annotation(self, value):
        for entry in self.identifiers():
            if entry.identifier == value:
                entry.remove_from_parent()
                return
        raise ValueError(f"No annotation found for '{value}'")

    def to_dict(self):
        return dict(self)

    def clear(self):
        for k in self:
            del self[k]
