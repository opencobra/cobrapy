"""Classes to handle standardized annotations and legacy dict-like annotations.

Standardized annotations enable users to associate MIRIAM-compliant annotations to model
components. They correspond to Controlled Vocabulary terms (CV terms), as described in
the SBML level 3 version 2 Core specification:
https://identifiers.org/combine.specifications:sbml.level-3.version-2.core.release-2
"""

import re
from collections import UserList
from collections.abc import Iterable as ABCIterable
from collections.abc import KeysView, MutableMapping
from typing import (  # TypeAlias, # Not supported in older python versions
    Any,
    Callable,
    Dict,
    FrozenSet,
    Generator,
    Iterable,
    Iterator,
    List,
    Optional,
    Pattern,
    Tuple,
    Union,
)

from cobra.core.metadata import Qualifier, Resource
from cobra.core.metadata.metadata import Metadata
from cobra.core.metadata.resource import QualifiersAlias, get_default_qualifier


# TypeAlias
StandardizedAnnotationInput = Union["StandardizedAnnotation", dict, str, Resource]


class StandardizedAnnotation:
    """Standardized annotation entry that defines a relation to MIRIAM resources.

    `StandardizedAnnotation` can be used to annotate cobra objects in a structural and
    standardized manner. This improves interoperability with other tools and neatly
    structures annotations. Therefore, standardized annotations should be preferred over
    custom annotations (`CustomAnnotation`), whenever possible.

    A `StandardizedAnnotation` object defines a set of MIRIAM (identifiers.org)
    resources and their relation to a modelling object. Relations are defined using
    the `Qualifier` enum, which defines a predefined set of qualifiers that either
    relate to the modelling object (e.g. `Qualifier.Modelling_is`) or relate to the
    biological object represented by the modelling object (e.g.
    `Qualifier.Biological_is`).

    For example, the cobra reaction representing a
    cytosolic transketolase reaction in Homo sapiens can be annotated with the
    resource "https://identifiers.org/reactome/R-HSA-163751" using the qualifier
    `Qualifier.Biological_is`, since the reactome entry R-HSA-163751 represents the
    biological reaction entity which is represented by the reaction object. Similarly,
    the RHEA resource "https://identifiers.org/rhea/RHEA:27628" can be added to the
    same `StandardizedAnnotation` object, since it relates to the same information about
    the cobra reaction and `Qualifier.Biological_is` is also applicable. The EC code
    ("https://identifiers.org/ec-code/2.2.1.1") for the same reaction should however be
    annotated using the `Biological_isVersionOf` qualifier, since the reaction can be
    seen as a version of the enzymatic activity represented by EC 2.2.1.1. Optionally, a
    standardized annotation object can have nested annotations. In the transketolase
    example, this could include a StandardizedAnnotation object with the qualifier
    `Qualifier.Biological_isDescribedBy` and the resource
    "https://identifiers.org/pubmed/9357955". This nested annotation adds additional
    information to the parent annotation, without changing its meaning.

    See `Biomodels Qualifiers
    <https://co.mbine.org/author/biomodels.net-qualifiers/>`_
    For a definition of all qualifiers, see `SBML Level 3, Version 2 Core, p 104
    <https://identifiers.org/combine.specifications:sbml.level-3.version-2.core.release-2>`_

    Parameters
    ----------
    resources: Resource, str, list or None, optional
        The resources as `Resource` objects or strings (URI format:
        https://identifiers.org/<namespace>/<id>). Default None.
    qualifier: Qualifier or str
        The qualifier for the relationship. Default `Qualifier.Biological_is`.
    annotations: list or None, optional
        List of StandardizedAnnotation objects that represent additional (nested)
        annotations of this object.
    """

    def __init__(
        self,
        resources: Optional[
            Union["Resource", str, Iterable[Union["Resource", str]]]
        ] = None,
        qualifier: Union[Qualifier, str] = Qualifier.Biological_is,
        annotations: Optional[Iterable["StandardizedAnnotation"]] = None,
    ):
        """Initialize a standardized annotation."""
        self._resources = []
        self.resources = resources
        self._qualifier: Qualifier = self.check_qualifier_type(qualifier)
        self._annotations: Optional[StandardizedAnnotationStore] = (
            StandardizedAnnotationStore.from_data(annotations)
            if annotations is not None
            else None
        )
        self._parent: Optional[StandardizedAnnotationStore] = None

    def _set_parent(self, parent):
        if self._parent is None or parent is None:
            self._parent = parent
        else:
            raise ValueError(
                "StandardizedAnnotation already has a parent. Create a new "
                "annotation if you would like to add an annotation to a second object."
            )

    def remove_from_parent(self) -> None:
        """Remove annotation from parent (`StandardizedAnnotationStore`).

        Raises
        ------
        ValueError
            If there is no known parent object.

        See Also
        --------
        StandardizedAnnotationStore.remove
        """
        if self._parent is None:
            raise ValueError(
                "Cannot remove annotation, since no object is associated with it."
            )
        self._parent.remove(self)

    def _remove_resource(self, resource):
        self.resources.remove(resource)

    @property
    def qualifier(self) -> Qualifier:
        """Get the qualifier.

        Returns
        -------
        Qualifier
        """
        return self._qualifier

    @qualifier.setter
    def qualifier(self, qualifier: Union[str, Qualifier]) -> None:
        """Set the qualifier.

        Parameters
        ----------
        qualifier: str or Qualifier

        See Also
        --------
        StandardizedAnnotation.check_qualifier_type()
        """
        self._qualifier = self.check_qualifier_type(qualifier)

    @property
    def resources(self) -> List["Resource"]:
        """Get the list of resources.

        Returns
        -------
        list of Resource objects
        """
        return self._resources

    @resources.setter
    def resources(self, resources: Iterable[Union[str, "Resource"]]) -> None:
        """Set the resources.

        Parameters
        ----------
        resources: list of str or Resource

        See Also
        --------
        StandardizedAnnotation.check_resource_type()
        """
        for idf in self._resources:
            idf._set_parent(None)
        self._resources = self.check_resource_type(resources)
        for idf in self._resources:
            idf._set_parent(self)

    def add_resources(self, resources: Iterable[Union[str, "Resource"]]) -> None:
        """Add resources to the standardized annotation.

        Parameters
        ----------
        resources: list of str or Resource objects
            List of resources to append to the existing ones.
        """
        resources = self.check_resource_type(resources)
        for resource in resources:
            resource._set_parent(self)
        self.resources.extend(resources)

    @property
    def uris(self) -> FrozenSet[str]:
        """Get the set of URIs represented by the resources of this annotation.

        Returns
        -------
        Set of URIs

        See Also
        --------
        all_uris
        Resource.uri
        """
        return frozenset({entry.uri for entry in self.resources})

    @property
    def all_uris(self) -> FrozenSet[str]:
        """Get all the URIs in this annotation, including nested annotations.

        Returns
        -------
        Set of URIs

        See Also
        --------
        uris
        Resource.uri
        """
        resources = {entry.uri for entry in self.resources}
        for entry in self.annotations:
            resources.update(entry.uris)
        return frozenset(resources)

    @property
    def annotations(self) -> "StandardizedAnnotationStore":
        """Get the nested annotations.

        Returns
        -------
        List of annotations
        """
        if self._annotations is None:
            self._annotations = StandardizedAnnotationStore()
        return self._annotations

    @annotations.setter
    def annotations(self, annotations: Iterable["StandardizedAnnotation"]) -> None:
        """Set the nested annotations.

        Parameters
        ----------
        annotations - list of StandardizedAnnotation objects.
        """
        self._annotations = StandardizedAnnotationStore.from_data(annotations)

    @staticmethod
    def check_resource_type(
        resources: Optional[
            Union[
                "Resource",
                str,
                Dict[str, str],
                Tuple[str, str],
                Iterable[Union["Resource", str, Dict[str, str], Tuple[str, str]]],
            ]
        ],
    ) -> List["Resource"]:
        """Check and parse resources.

        Parameters
        ----------
        resources: Resource, str, dict, tuple, or list thereof, optional
            Input data to check if it is or can be transformed to Resource objects.
            String must start with http:// or https:// to be acceptable.
            Dictionary must match the format required by from_data.
            If None is given, an empty StandardizedAnnotation object is returned.
            No parsing of identifiers/URIs is done, perhaps in future versions.

        Returns
        -------
        list of Resource objects

        See Also
        --------
        Resource.from_dict
        """
        if resources is None:
            return []
        elif isinstance(resources, (Resource, str, dict)):
            return [Resource.from_data(resources, strict=False)]
        elif isinstance(resources, ABCIterable):
            return [
                Resource.from_data(resource, strict=False) for resource in resources
            ]
        else:
            raise TypeError(
                f"Allowed types for resources are Resource, str, or a list thereof,"
                f"not {type(resources)}: {resources}"
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
        if isinstance(qual, str) and qual not in Qualifier._value2member_map_:
            raise TypeError(f"{qual} is not a supported enum Qualifier")
        elif isinstance(qual, Qualifier):
            return qual
        elif isinstance(qual, str):
            return Qualifier._value2member_map_[qual]
        else:
            raise TypeError(
                f"Allowed types for StandardAnnotation qualifiers are Qualifier or"
                f"str member of the Qualifier enum {type(qual)}, {qual}"
            )

    def to_dict(self) -> Dict:
        """Convert annotation to a python dict.

        Returns
        -------
        dict:
            A dict that has up to three keys
            "qualifier" - the qualifier as a string
            "resources" - the resources as list
            "annotations" (optionally) - the nested annotations as a list
        """
        d = {
            "qualifier": self.qualifier.value,
            "resources": [resource.uri for resource in self.resources],
        }
        if self.annotations:
            d["annotations"] = self.annotations.to_list_of_dicts()
        return d

    def to_tuples(self) -> List[Tuple[str, str]]:
        """Convert the annotation to a list of tuples of namespace-identifier pairs.

        This does not contain the qualifier or nested annotations.

        Returns
        -------
        List of namespace-identifier pairs as tuples
        """
        return [
            (resource.namespace, resource.identifier)
            for resource in self.resources
            if resource.namespace is not None and resource.identifier is not None
        ]

    def to_records(self) -> List[Dict[str, Any]]:
        """Convert the annotation to a list of dictionaries that represent resources.

        Each entry in the list represents a resource associated with this object,
        either directly or as nested annotation. The resulting list is thus a flattened
        representation of all resources. The hierarchical information is included
        using the "annotation_group" and "parent_group" entries in the dictionaries.
        Each annotation group represents a single StandardizedAnnotation object (without
        its nested annotations). If a resource is found in a nested annotation, the
        "parent_group" value will be set to the "annotation_group" value of its parent
        StandardizedAnnotation object.

        Returns
        -------
        list of dict objects
            Each dict represents a Resource records.

        See Also
        --------
        StandardizedAnnotationStore.to_records
        """
        records, _ = self._to_records()
        return records

    def _to_records(
        self, group_counter: int = 1, parent_group: int = 0
    ) -> Tuple[List[Dict[str, Any]], int]:
        records = []
        for resource in self.resources:
            entry = {
                "qualifier": self.qualifier.value,
                "uri": resource.uri,
                "namespace": resource.namespace,
                "identifier": resource.identifier,
                "annotation_group": group_counter,
                "parent_group": parent_group,
            }
            records.append(entry)
        if self.annotations:
            new_records, group_counter = self.annotations._to_records(
                group_counter=(group_counter + 1), parent_group=group_counter
            )
            records.extend(new_records)
        else:
            group_counter = group_counter + 1
        return records, group_counter

    @classmethod
    def from_dict(cls, data_dict: Dict) -> "StandardizedAnnotation":
        """Generate a StandardizedAnnotation object from a python dict.

        Parameters
        ----------
        data_dict: dict
            A dict that has up to three keys
            "qualifier" - the qualifier as a string, optional. If not present, the
            qualifier is set to `Qualifier.Biological_is`.
            "resources" - the resources as a list of strings, optional.
            "annotations" - the nested annotations as a list of dicts, optional.

        Returns
        -------
        StandardizedAnnotation

        See Also
        --------
        to_dict

        """
        return cls(
            resources=data_dict.get("resources", None),
            qualifier=data_dict.get("qualifier", Qualifier.Biological_is),
            annotations=data_dict.get("annotations", None),
        )

    def __eq__(self, other: Any) -> bool:
        """Compare StandardizedAnnotation to another object and determine equality.

        If a dict is given, it is transformed to `StandardizedAnnotation` before
        comparison. Will return False for any other type. The order of the resources
        is ignored.

        Parameters
        ----------
        other

        Returns
        -------
        bool
            False if other is not StandardizedAnnotation or dict.
            False if qualifiers, resources or nested annotations are different.
            True otherwise.

        See Also
        --------
        StandardizedAnnotation.from_dict()
        """

        if isinstance(other, dict):
            return self == StandardizedAnnotation.from_dict(other)
        if isinstance(other, StandardizedAnnotation):
            if self.qualifier != other.qualifier:
                return False
            if len(self.resources) != len(other.resources):
                return False
            for idf in self.resources:
                if idf not in other.resources:
                    return False
            if self.annotations != other.annotations:
                return False
            return True
        return False

    def __repr__(self) -> str:
        """Return the StandardizedAnnotation as str with module and class.

        Returns
        -------
        str
        """
        return (
            f"{self.__class__.__module__}.{self.__class__.__qualname__}"
            f"({self.to_dict()})"
        )

    def _repr_html_(self) -> str:
        """Return the StandardizedAnnotation as HTML string.

        Returns
        -------
        str
            HTML formatted string
        """
        # TODO: Fix this HTML
        cols = {
            "uri": "<td><strong>URI</strong></td>",
            "namespace": "<td><strong>Namespace</strong></td>",
            "identifier": "<td><strong>Identifier</strong></td>",
            "address": "<td><strong>Memory address</strong></td>",
        }
        n = len(self.resources)
        s = f"""<table>
            <tr>
                <td><strong>Qualifier</strong></td>
                <td colspan='{n}'
                    style='border-bottom: solid black 1px; text-align: center;'>
                    {self.qualifier}
                </td>
            </tr>"""
        for resource in self.resources:
            cols["uri"] += f"<td>{resource.uri}</td>"
            cols["namespace"] += f"<td>{resource.namespace}</td>"
            cols["identifier"] += f"<td>{resource.identifier}</td>"
            cols["address"] += f"<td>{id(resource):#x}</td>"
        for v in cols.values():
            s += f"<tr>{v}</tr>"
        if self.annotations:
            s += f"""<tr>
                <td><strong>Nested<br />annotations</strong></td>
                <td colspan='{n}'>"""
            for annotation in self.annotations:
                s += annotation._repr_html_()
            s += "</td></tr>"

        s += "</table>"
        return s

    def __deepcopy__(self, memo: Optional[dict] = None):
        """Copy the StandardizedAnnotation efficiently with memo.

        Parameters
        ----------
        memo: dict
            Automatically passed parameter, dict of already copied items.

        Returns
        -------
        StandardizedAnnotation
            A new annotation that is a deep copy of the original.
        """

        new = StandardizedAnnotation(
            qualifier=self.qualifier,
        )
        # memo[id(self)] = new
        new.resources = [r.__deepcopy__() for r in self._resources]
        if self._annotations is not None:
            new._annotations = self._annotations.__deepcopy__()
        return new

    def copy(self) -> "StandardizedAnnotation":
        """Copy the annotation and all its nested annotations.

        Returns
        -------
        StandardizedAnnotation
            A new annotation that is a deep copy of the original.
        """

        return self.__deepcopy__()


class StandardizedAnnotationStore(UserList):
    """A list-like object that stores a collection of `StandardizedAnnotation` objects.

    Stores a collection of standardized annotations that define the relation  of
    MIRIAM-compliant resources to a cobra modelling object. In practice, this class is
    automatically instantiated as the `standardized` attribute of the `Metadata` class,
    which can be accessed through `object.metadata.standardized`. In addition, nested
    annotations in a `StandardizedAnnotation` object also make use of the
    `StandardizedAnnotationStore` class.

    Parameters
    ----------
    data: list of StandardizedAnnotation, dict, str or Resource objects
        List of standardized annotations to initialize the store with. Dictionaries
        are converted to standardized annotations using
        `StandardizedAnnotation.from_dict`. Strings are interpreted as identifier
        URIs and together with other Resource objects stored in a new
        `StandardizedAnnotation` with `Qualifier.Biological_is` as qualifier.
    """

    def __init__(
        self,
        data: Optional[
            Iterable[Union[StandardizedAnnotation, Dict, Resource, str]]
        ] = None,
    ):
        """Initialize a standardized annotation store."""
        self._take_ownership_of_resources: bool = getattr(
            self, "_take_ownership_of_resources", True
        )
        if data is None:
            data = []

        checked_data = [
            filtered_entry
            for entry in data
            if (not isinstance(entry, (str, Resource)))
            and (filtered_entry := self._check_standardized_annotation(entry))
            is not None
        ]
        # str and Resource instances are handled separately and added as one
        # Standardized annotation instance with the default qualifier (Biological_is).
        if no_qualifier_data := [
            entry for entry in data if isinstance(entry, (str, Resource))
        ]:
            checked_data.insert(0, StandardizedAnnotation(resources=no_qualifier_data))
        if self._take_ownership_of_resources:
            for entry in checked_data:
                entry._set_parent(self)
        super().__init__(checked_data)

    @staticmethod
    def _check_standardized_annotation(
        ann: Optional[Union[StandardizedAnnotation, Dict, str, Resource]],
    ) -> Optional["StandardizedAnnotation"]:
        if ann is None:
            return None
        if isinstance(ann, StandardizedAnnotation):
            return ann
        elif isinstance(ann, (str, Resource)):
            return StandardizedAnnotation(ann)
        elif isinstance(ann, dict):
            return StandardizedAnnotation.from_dict(ann)
        else:
            raise TypeError(
                f"Allowed types for StandardizedAnnotation are str and"
                f"StandardizedAnnotation, not {type(ann)}: {ann}"
            )

    @staticmethod
    def from_data(
        data: Optional[
            Union[
                Iterable[Union[str, Dict, "StandardizedAnnotation"]],
                str,
                Dict,
                "StandardizedAnnotation",
                "StandardizedAnnotationStore",
            ]
        ],
    ) -> "StandardizedAnnotationStore":
        """Create a StandardizedAnnotationStore object from given data.

        Parameters
        ----------
        data: StandardizedAnnotation, dict, str or Resource, or list thereof, or
        StandardizedAnnotationStore or None
            A standardized annotation or a list of standardized annotations to use to
            create a store with. Dictionaries are converted to standardized annotations
            using `StandardizedAnnotation.from_dict`. Strings are interpreted as
            identifier URIs and together with other Resource objects stored in a new
            `StandardizedAnnotation` with `Qualifier.Biological_is` as qualifier.
            If data is already a StandardizedAnnotationStore, this object will simply be
            returned.

        Returns
        -------
        StandardizedAnnotationStore

        Raises
        ------
        TypeError
        """
        if data is None:
            return StandardizedAnnotationStore()
        elif isinstance(data, StandardizedAnnotationStore):
            return data
        elif isinstance(data, (StandardizedAnnotation, dict, str)):
            if not data:
                return StandardizedAnnotationStore()
            return StandardizedAnnotationStore([data])
        elif isinstance(data, ABCIterable):
            return StandardizedAnnotationStore(data)
        else:
            raise TypeError(f"Invalid format for StandardizedAnnotationStore: '{data}'")

    def to_list_of_dicts(self) -> List[dict]:
        """Convert the StandardizedAnnotationStore to a list of python dicts.

        Returns:
        -------
        list:
            a list where each item is a dict representing a standardized annotation
            object, as created by StandardizedAnnotation.to_dict(). Mainly used for JSON
            and YAML export.

        See Also
        --------
        StandardizedAnnotation.to_dict
        """
        return [annotation.to_dict() for annotation in self]

    def to_records(self):
        """Convert the store to a list of dictionaries that represent resources.

        Each entry in the list represents an resources associated with one of the
        standardized annotation in this object, either directly or as nested annotation.
        The resulting list is thus a flattened representation of all resources.
        The hierarchical information is included using the "annotation_group" and
        "parent_group" entries in the dictionaries. Each annotation group represents a
        single StandardizedAnnotation object (without its nested annotations). If a
        resource is found in a nested annotation, the "parent_group" value will be set
        to the "annotation_group" value of its parent StandardizedAnnotation object.

        Returns
        -------
        List of dictionaries representing Resource records

        See Also
        --------
        StandardizedAnnotation.to_records
        """
        records, _ = self._to_records()
        return records

    def _to_records(
        self, group_counter: int = 1, parent_group: int = 0
    ) -> Tuple[List[Dict], int]:
        records = []
        for entry in self.data:
            new_l, group_counter = entry._to_records(
                group_counter=group_counter, parent_group=parent_group
            )
            records.extend(new_l)
        return records, group_counter

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
            entry = StandardizedAnnotation(qualifier=qualifier, resources=[])
            if self._take_ownership_of_resources:
                entry._set_parent(self)
            self.data.insert(0, entry)
        return entry

    def add(
        self,
        annotations: Union[
            "StandardizedAnnotationStore",
            StandardizedAnnotation,
            Dict,
            str,
            Resource,
            List[Union[StandardizedAnnotation, Dict, str, Resource]],
        ],
    ) -> None:
        """Add one or multiple standardized annotations to the store.

        Dictionaries are converted to standardized annotations using
        `StandardizedAnnotation.from_dict`. Strings are interpreted as identifier
        URIs and together with other Resource objects stored in a new
        `StandardizedAnnotation` with `Qualifier.Biological_is` as qualifier.

        If a list is passed, this method is equivalent to
        `StandardizedAnnotationStore.extend` and if a single annotation is passed, this
        method is equivalent to `StandardizedAnnotationStore.append`.

        Parameters
        ----------
        annotations : StandardizedAnnotation, dict, str or Resource, or list thereof
            Single or multiple annotations to add to the store.

        See Also
        --------
        append
        extend
        """
        if isinstance(annotations, (StandardizedAnnotation, dict, str, Resource)):
            self.append(annotations)
        else:
            self.extend(annotations)

    def append(self, item: Union[StandardizedAnnotation, Dict, str, Resource]) -> None:
        """Append a single annotation to the end of the store.

        Parameters
        ----------
        item: StandardizedAnnotation, dict, str or Resource
            Annotation to append to the standardized annotation store.
        """
        self.extend([item])

    def extend(
        self,
        other: Union[
            "StandardizedAnnotationStore",
            Iterable[Union[StandardizedAnnotation, Dict, str, Resource]],
        ],
    ) -> None:
        """Extend store by appending elements from the iterable.

        Parameters
        ----------
        iterable : Iterable
            Annotations to add to the store.
        """
        if isinstance(other, StandardizedAnnotationStore):
            self.extend(other.data)
        elif isinstance(other, Iterable):
            checked_data = [
                checked_item
                for item in other
                if (checked_item := self._check_standardized_annotation(item))
                is not None
            ]
            if self._take_ownership_of_resources:
                for d in checked_data:
                    d._set_parent(self)
            self.data.extend(checked_data)

    def remove(self, item: "StandardizedAnnotation"):
        """Remove a standardized annotation from the store.

        Parameters
        ----------
        item: StandardizedAnnotation
            Annotation to remove from the store.
        """
        self.data.remove(item)
        item._set_parent(None)

    @property
    def resources(self) -> FrozenSet[Resource]:
        """Get resources.

        Returns
        -------
        FrozenSet
            a set of resources in the standardized annotation store, not including
            resources of nested annotations.
        """
        resources = set()
        for entry in self.data:
            resources.update(entry.resources)
        return frozenset(resources)

    @property
    def all_resources(self) -> FrozenSet[Resource]:
        """Get all resources, including resources in nested annotations.

        Returns
        -------
        FrozenSet
            a set of all resources in the standardized annotation store, including
            resources of nested annotations.
        """
        resources = set()
        for entry in self.data:
            resources.update(entry.resources)
            if entry.annotations:
                resources.update(entry.annotations.all_resources)
        return frozenset(resources)

    def resources_for(
        self,
        namespace: Optional[Union[str, List[str]]] = None,
        qualifier: Optional[
            Union[
                Qualifier,
                QualifiersAlias,
                List[Union[Qualifier, QualifiersAlias]],
            ]
        ] = None,
        nested: bool = False,
    ) -> List[Resource]:
        """Get resources for a namespace or qualifier.

        Filter annotations based on their qualifier and its resources on their
        namespace. Optionally also return and filter nested resources on their
        namespace.

        Parameters
        ----------
        namespace: None, str or list of str, optional
            One or multiple namespaces to filter resources with. Selects a resource when
            its namespace matches any of the provided namespaces. If it is set to None,
            no filtering based on the namespace will be performed. Default None.
        qualifier: None, Qualifier, QualifiersAlias or list of Qualifier/QualifiersAlias
            One or multiple qualifiers to filter annotations with. Selects an annotation
            when its qualifier matches any of the provided qualifiers. If it is set to
            None, no filtering based on qualifiers will be performed. Nested annotations
            are never filtered based on their qualifier. Default None.
        nested: bool
            Whether to return resources from nested annotations. Nested annotations are
            selected when `nested` is True and the top-level annotation is selected
            based on its qualifier. I.e. nested annotations are not filtered on their
            own qualifier. Conversely, resources in nested annotations are selected
            based on their namespace. Default False.

        Returns
        -------
        list of Resource objects
        """
        # TODO: Examples
        namespace_sel, qualifier_sel = True, True
        if namespace is None:
            namespace_sel = False
            namespace = []
        if not isinstance(namespace, list):
            namespace = [namespace]

        if qualifier is None:
            qualifier_sel = False
            qualifier = []
        if not isinstance(qualifier, list):
            qualifier = [qualifier]
        qualifiers_and_aliases = qualifier
        qualifier = []
        for q in qualifiers_and_aliases:
            if isinstance(q, Qualifier):
                qualifier.append(q)
            elif isinstance(q, QualifiersAlias):
                qualifier.extend(q)
            else:
                raise TypeError(
                    "Qualifiers should have type Qualifier or QualifiersAlias,"
                    f"not {type(q)}"
                )

        qualifier_set = set(qualifier)
        namespace_set = set(namespace)
        resources = []
        for annotation in self:
            if qualifier_sel and annotation.qualifier not in qualifier_set:
                continue
            if not namespace_sel:
                resources.extend(annotation.resources)
            else:
                for resource in annotation.resources:
                    if resource.namespace is None:
                        continue
                    if resource.namespace in namespace_set:
                        resources.append(resource)
            if nested and annotation.annotations:
                # Do not select for qualifiers in nested annotations
                resources.extend(
                    annotation.annotations.resources_for(
                        namespace=list(namespace_set) if namespace_sel else None,
                        nested=nested,
                    )
                )
        return resources

    @property
    def uris(self) -> FrozenSet[str]:
        """Get URIs.

        Returns
        -------
        FrozenSet
            A set of URIs in the standardized annotation store, not including URIs of
            nested annotations.
        """
        resources = set()
        for entry in self.data:
            resources.update(entry.uris)
        return frozenset(resources)

    @property
    def all_uris(self) -> FrozenSet[str]:
        """Get all URIs, including URIs of nested annotations.

        Returns
        -------
        FrozenSet
            A set of URIs in the standardized annotation store, including URIs of
            nested annotations.
        """
        resources = set()
        for entry in self.data:
            resources.update(entry.uris)
            if entry.annotations:
                resources.update(entry.annotations.all_uris)
        return frozenset(resources)

    @property
    def qualifiers(self) -> FrozenSet[Qualifier]:
        """Get qualifiers of annotations in the store, not including nested annotations.

        Returns
        -------
        FrozenSet
            A set of qualifiers in the standardized annotation store, not including
            qualifiers of nested annotations.
        """
        qualifier_set = set()
        for entry in self.data:
            qualifier_set.add(entry.qualifier)
        return frozenset(qualifier_set)

    @property
    def all_qualifiers(self) -> FrozenSet[Qualifier]:
        """Get all qualifiers of annotations in the store, including nested annotations.

        Returns
        -------
        FrozenSet
            A set of qualifiers in the standardized annotation store, including
            qualifiers of nested annotations.
        """
        qualifier_set = set()
        for entry in self.data:
            qualifier_set.add(entry.qualifier)
            if entry.annotations:
                qualifier_set.update(entry.annotations.all_qualifiers)
        return frozenset(qualifier_set)

    def __iter__(self) -> Iterator[StandardizedAnnotation]:
        """Get an iterator for the standardized annotations in the store.

        Returns
        -------
        Iterator
        """
        return iter(self.data)

    def __len__(self) -> int:
        """Get the number of standardized annotations in the store.

        Returns
        -------
        int
        """
        return len(self.data)

    def query(
        self,
        search_function: Union[str, Pattern, Callable],
        attribute: Union[str, None] = None,
    ) -> "StandardizedAnnotationList":
        """Query the annotation store for matching StandardizedAnnotation objects.

        Parameters
        ----------
        search_function : a string, regular expression or function
            Used to find the matching elements in the store.
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
        StandardizedAnnotationList
            A list-like collection of StandardizedAnnotation objects which match the
            query.

        Examples
        --------
        >>> from cobra.io import load_model
        >>> model = load_model('iJO1366')
        >>> model.metadata.standardized.query('Biological_', 'qualifier')
        >>> import re
        >>> regex = re.compile('^Modelling')
        >>> model.annotation.standardized.query(regex, 'qualifier')
        """

        # TODO: Clean up this whole method.
        def select_attribute(
            x: StandardizedAnnotation,
        ) -> Union[StandardizedAnnotation, Resource, Qualifier, set]:
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
                    annotation
                    for annotation in self
                    if (
                        regex_searcher.findall(select_attribute(annotation).name) != []
                        or regex_searcher.findall(select_attribute(annotation).value)
                        != []
                    )
                ]
            elif attribute == "resources":
                matches = [
                    annotation
                    for annotation in self.data
                    if any(
                        regex_searcher.findall(res.uri)
                        for res in select_attribute(annotation)
                    )
                ]
            else:
                matches = [
                    annotation
                    for annotation in self.data
                    if regex_searcher.findall(annotation.qualifier.name) != []
                    or regex_searcher.findall(annotation.qualifier.value) != []
                    or any(
                        regex_searcher.findall(res.uri) for res in annotation.resources
                    )
                ]
        except TypeError as err:
            print(err)
            matches = [
                annotation
                for annotation in self.data
                if search_function(select_attribute(annotation))
            ]

        return StandardizedAnnotationList(matches)

    def __setitem__(self, key: int, value: StandardizedAnnotationInput) -> None:
        """Set item in the store at the provided index.

        Removes the current annotation at the provided index and replaces it with the
        provided annotation.

        Parameters
        ----------
        key: int
        value: StandardizedAnnotation, dict, str or Resource
        """
        checked_value = self._check_standardized_annotation(value)
        if checked_value is None:
            raise TypeError("Value cannot be None.")
            # TODO: Elaborate (or automatically delete when None)
        if self._take_ownership_of_resources:
            self.data[key]._set_parent(None)
            checked_value._set_parent(self)
        UserList.__setitem__(self, key, checked_value)

    def __getitem__(
        self,
        key: Union[
            int,
            Qualifier,
            QualifiersAlias,
            List[Union[int, Qualifier, QualifiersAlias]],
        ],
    ) -> Union[StandardizedAnnotation, "StandardizedAnnotationList"]:
        """Access standardized annotations by integer index or qualifier.

        Parameters
        ----------
        key: int, Qualifier, QualifiersAlias or list thereof
            If `key` is an integer, the `StandardizedAnnotation` at that position in the
            list will be returned. When a list of integers is provided, a
            `StandardizedAnnotationList` with the corresponding annotations is returned.
            When one or more Qualifier or `QualifiersAlias` enums are provided, a
            `StandardizedAnnotationList` of all annotations (not nested) with any of
            those qualifiers is returned.
        """
        if isinstance(key, int):
            return self.data[key]
        if isinstance(key, (Qualifier, QualifiersAlias)):
            key = [key]
        if not isinstance(key, list):
            raise TypeError(f"Indexed using key of wrong type: {type(key)}")
        selection = []
        for k in key:
            if isinstance(k, int):
                selection.append(self.data[k])
            elif isinstance(k, Qualifier):
                selection.extend(
                    annotation for annotation in self if annotation.qualifier == k
                )
            elif isinstance(k, QualifiersAlias):
                selection.extend(
                    annotation for annotation in self if annotation.qualifier in k
                )
            else:
                raise TypeError(f"Indexed using a key of wrong type: {type(key)}")
        return StandardizedAnnotationList(selection)

    def __eq__(self, other: Any) -> bool:
        """Compare two standardized annotation stores and determine equality.

        Equality is defined as them having the same data, but not necessarily being the
        same object. If the given item is not a StandardizedAnnotationStore, list or
        dict, this function will return False.

        Parameters
        ----------
        other

        Returns
        -------
        bool: True if the data matches, False otherwise
        """
        if isinstance(other, (ABCIterable, dict)) and not isinstance(
            other, StandardizedAnnotationStore
        ):
            return self.__eq__(StandardizedAnnotationStore.from_data(other))
        if not isinstance(other, StandardizedAnnotationStore):
            return False
        if len(self.data) != len(other.data):
            return False
        for other_entry in other.data:
            if other_entry not in self.data:
                return False
        return True

    def _repr_html_(self) -> str:
        """Convert StandardizedAnnotationStore to HTML.

        Returns
        -------
        str
            HTML representation of the annotation store.
        """
        entries = [annotation._repr_html_() for annotation in self]
        return f"""StandardizedAnnotationStore{"<p>".join(entries)}"""

    def __deepcopy__(self, memo: Optional[dict] = None):
        """Copy the StandardizedAnnotationStore efficiently with memo.

        Parameters
        ----------
        memo: dict
            Automatically passed parameter, dict of already copied items.

        Returns
        -------
        StandardizedAnnotationStore
            A new annotation store that is a deep copy of the original.
        """

        new = StandardizedAnnotationStore()
        # memo[id(self)] = new
        new.extend([ann.__deepcopy__() for ann in self])
        return new

    def copy(self) -> "StandardizedAnnotationStore":
        """Copy the annotation store and all its annotations.

        Returns
        -------
        StandardizedAnnotationStore
            A new store that is a deep copy of the original store.
        """

        return self.__deepcopy__()


class StandardizedAnnotationList(StandardizedAnnotationStore):
    """Class to create lists of StandardizedAnnotation objects.

    This class is very similar to the StandardizedAnnotationStore class, which stores
    all standardized annotation objects of a cobrapy object. The difference is that this
    class does not take ownership (becomes parent of) its resources. It can therefore be
    used to create lists of StandardizedAnnotation objects of different cobrapy objects,
    or create views of a subset of the annotations of a single object.

    Parameters
    ----------
    data: list of StandardizedAnnotation, dict, str or Resource objects
        List of standardized annotations to initialize the store with. Dictionaries
        are converted to standardized annotations using
        `StandardizedAnnotation.from_dict`. Strings are interpreted as identifier
        URIs and together with other Resource objects stored in a new
        `StandardizedAnnotation` with `Qualifier.Biological_is` as qualifier.

    """

    def __init__(
        self,
        data: Optional[
            Iterable[Union[StandardizedAnnotation, Dict, Resource, str]]
        ] = None,
    ):
        """Initialize a standardized annotation store."""

        self._take_ownership_of_resources: bool = False
        super(StandardizedAnnotationList, self).__init__(data)


class SimplifiedAnnotationInterface(MutableMapping):
    """Class to interface with metadata using a dict-like interface.

    This class is used to maintain compatibility with older cobrapy versions. It is
    typically accessed through an objects annotation attribute. It allows a user to get
    and set standardized annotations through a dict-like interface. When reading
    existing annotations, qualifiers are ignored and resources are pooled. When
    setting new annotations, qualifiers are set based on defaults (typically
    `Qualifiers.Biological_is`).
    This class is automatically instantiated as the `annotation` attribute of cobrapy
    objects.

    Warnings
    --------
    * This is not the preferred method to access annotations, since information and
        hierarchy is lost in this interface.
    * This interface was added to not break existing code, for new code
        `object.metadata.standardized` should be preferred.
    * Editing existing annotations using this interface can cause the annotations to
        become less organized.

    Parameters
    ----------
    metadata: Metadata
        Metadata object where annotations will be stored and retrieved from.
    """

    def __init__(self, metadata: Metadata) -> None:
        """Initialize the simplified annotation interface using a `Metadata` object."""
        self._metadata: Metadata = metadata

    def add(
        self,
        data: Optional[
            Union[
                Dict,
                str,
                "SimplifiedAnnotationInterface",
                Tuple[str, Union[str, List[str]]],
                Resource,
                List[Union[str, Resource, Tuple[str, Union[str, List[str]]]]],
            ]
        ] = None,
    ) -> None:
        """Add an annotation.

        Parameters
        ----------
        data: str, tuple, Resource or list thereof, or dict or
        SimplifiedAnnotationInterface
            Add resources as annotations, using default qualifiers. Tuples are
            interpreted as namespace-identifiers pairs, strings should be valid URIs and
            if a dictionary is provided, its keys should represent namespaces and its
            values identifiers.

        Raises
        ------
        TypeError
            If values could not be converted to Resource objects.
        ValueError
            If tuples of lengths other than 2 were provided.
        """
        if data is None:
            data = []

        if isinstance(data, SimplifiedAnnotationInterface):
            data = data.to_dict()

        if isinstance(data, dict):
            data = list(data.items())

        if isinstance(data, (str, tuple)):
            data = [data]

        if not isinstance(data, list):
            raise TypeError(
                "The supplied annotations were not of type List, or could "
                "not be converted to a list."
            )

        annotations = {}
        for entry in data:
            if isinstance(entry, tuple):
                if len(entry) != 2:
                    raise ValueError(
                        "Only tuples of length 2 can be converted to annotations."
                    )
                expanded_entries = []
                if isinstance(entry[1], list):
                    expanded_entries.extend([(entry[0], v) for v in entry[1]])
                else:
                    expanded_entries.append(entry)
                entry = expanded_entries
            elif isinstance(entry, str):
                entry = [entry]
            else:
                raise TypeError("Entry could could not be converted to an Resource.")
            for x in entry:
                if isinstance(x, tuple) and x[0].lower() == "sbo":
                    self._metadata.sbo = x[1]
                    continue
                if not isinstance(x, Resource):
                    x = Resource.from_data(x, strict=False)
                if x.namespace is None:
                    raise ValueError(f"Could not determine namespace of resource {x}.")

                qualifier = get_default_qualifier(x.namespace)
                if qualifier not in annotations:
                    annotations[qualifier] = []
                annotations[qualifier].append(x)

        for qualifier, resources in annotations.items():
            ann = self._metadata.standardized._find_first_or_create_by_qualifier(
                qualifier
            )
            ann.add_resources(resources)

    @property
    def sbo(self) -> str:
        """Get the SBO term of the annotations.

        Returns
        -------
        str

        See Also
        --------
        Metadata.sbo
        """
        return self._metadata.sbo

    @sbo.setter
    def sbo(self, value: str) -> None:
        """Set the SBO term of the annotations.

        Parameters
        ----------
        value: str

        See Also
        --------
        Metadata.sbo
        """
        self._metadata.sbo = value

    def __getitem__(self, key: str) -> List[str]:
        """Get resources for a given namespace.

        Collects all standardized annotations for namespace `key` and returns a list of
        all resources as strings, or raise IndexError if none were found.

        Parameters
        ----------
        key: str
            Namespace of the resources.

        Returns
        -------
        list of str
            List of resources as strings.

        Raises
        ------
        IndexError
            If no resources were found for the given namespace.
        """
        if not isinstance(key, str):
            raise TypeError("Index should be of type str.")
        key = key.lower()
        if key == "sbo":
            return [self._metadata.sbo]

        results = []
        for ann in self._metadata.standardized:
            for resource in ann.resources:
                if (v := resource.identifier) is not None and resource.namespace == key:
                    results.append(v)

        if results:
            # Deduplicate and sort results to have consistent results and make
            # comparisons easier. E.g. __eq__(...) relies on this.
            return list(sorted(set(results)))
        else:
            raise IndexError(f"No resources found for namespace '{key}'")

    def get(self, key: str, default: Any = None):
        """Get resources for a given namespace.

        Performs indexing in the same manner as `__get__`, except it will not raise
        IndexError when the index does not exist, but rather returns a default value.

        Parameters
        ----------
        key: str
            Namespace of the resources.
        default: optional
            Default value to return when index is not found. Default None.

        Returns
        -------
        list of str or `default`
            List of resources as strings or the `default` value.
        """
        # Example of usage:
        # https://github.com/opencobra/memote/blob/develop/src/memote/support/thermodynamics.py
        try:
            return self[key]
        except IndexError:
            return default

    def __delitem__(self, key: str) -> None:
        """Delete all resources for a given namespace.

        Deletes all resources for namespace `key` from their `StandardizedAnnotation`
        objects or raise IndexError if none were found.

        Parameters
        ----------
        key: str
            Namespace of the resources.

        Raises
        ------
        IndexError
            If no resources were found for the given namespace.
        """
        if not isinstance(key, str):
            raise TypeError("Index should be of type str.")

        key = key.lower()

        if key == "sbo":
            self._metadata.sbo = ""
            return

        deleted_any = False
        for ann in self._metadata.standardized:
            for resource in list(ann.resources):
                if resource.namespace is not None and resource.namespace == key:
                    resource.remove_from_parent()
                    deleted_any = True
        if not deleted_any:
            raise IndexError(f"Could not find annotations for f'{key}')")

    def __setitem__(self, key: str, value: Union[str, List[str]]) -> None:
        """Set resources for a given namespace.

        Removes all existing resources with namespace `key` and inserts
        the provided resources.

        Parameters
        ----------
        key: str
            Namespace of the resources.
        value: str or list of str
            Resources to set for the provided namespace.
        """
        if not isinstance(key, str):
            raise TypeError("Index should be of type str.")
        key = key.lower()

        try:
            del self[key]
        except IndexError:
            pass

        self.add({key: value})

    def __eq__(self, other: Any) -> bool:
        """Determine equality between the simplified annotations and another object.

        If the other object is of a type other than dict or
        `StandardizedAnnotationInterface`, the objects are considered not equal.
        Otherwise, objects are equal if there dictionary representation is equal.

        Parameters
        ----------
        other

        Returns
        -------
        bool
        """
        self_dict = self.to_dict()
        if isinstance(other, dict):
            other_dict = other
        elif isinstance(other, SimplifiedAnnotationInterface):
            other_dict = other.to_dict()
        else:
            return False
        return self_dict == other_dict

    def objects(self) -> Generator[Tuple[str, List[Resource]], None, None]:
        """Get a generator for all pairs of namespace and Resource objects.

        Creates a generator that can be used to iterate over the data as pairs (tuples)
        of a namespace and a list of its corresponding Resource objects. This method
        is very similar to the `items` method, except that it yields Resource objects
        instead of strings.

        Returns
        -------
        generator
            Yields pairs of namespace and a list of Resource objects.

        See Also
        --------
        items
        """
        visited_namespaces = set("sbo")
        while True:
            current_namespace = None
            resources = []
            for entry in self._metadata.standardized:
                for resource in entry.resources:
                    if resource.namespace is None:
                        continue
                    if resource.namespace in visited_namespaces:
                        continue
                    if current_namespace is None:
                        current_namespace = resource.namespace
                        resources.append(resource)
                    else:
                        if resource.namespace == current_namespace:
                            resources.append(resource)
            if current_namespace is not None:
                visited_namespaces.add(current_namespace)
            resources = list(sorted(set(resources), key=lambda x: x.identifier))
            if resources:
                yield (current_namespace, resources)
            else:
                break
        if sbo_term := self._metadata.sbo:
            yield ("sbo", [Resource.from_data(("sbo", sbo_term))])

    def items(self) -> Iterator[Tuple[str, List[str]]]:
        """Get a generator for all pairs of namespace and identifiers.

        Creates a generator that can be used to iterate over the data as pairs (tuples)
        of a namespace and a list of its corresponding identifiers. This method is very
        similar to the `objects` method, except that it yields strings instead of
        Resource objects.

        Returns
        -------
        generator
            Yields pairs of namespace and a list of identifiers.

        See Also
        --------
        objects
        """
        for k, v in self.objects():
            yield (k, [x.identifier for x in v])

    def __iter__(self):
        """Get an iterator for the namespaces (keys) of the annotations.

        Returns
        -------
        Iterator
        """

        for k, _ in self.objects():
            yield k

    def keys(self):
        """Get a view of the namespaces (keys) of the annotations.

        Returns
        -------
        Iterator
        """

        return KeysView(self)

    def values(self):
        """Get a generator for the values of the annotations.

        The values are lists of identifiers with the same namespace.

        Returns
        -------
        Generator

        See Also
        --------
        identifiers
        tuples
        """
        for _, v in self.items():
            yield v

    def identifiers(self):
        """Get a generator for the individual identifiers of the annotations.

        Returns
        -------
        Generator

        See Also
        --------
        values
        tuples
        """
        for _, v in self.items():
            for identifier in v:
                yield identifier

    def tuples(self):
        """Get a generator for the individual identifiers as namespace-identifier pairs.

        Returns
        -------
        Generator

        See Also
        --------
        values
        identifiers
        """
        for k, v in self.items():
            for identifier in v:
                yield (k, identifier)

    @property
    def number_of_resources(self) -> int:
        """The number of individual resources.

        Is different from the length of the object, since a single namespace can have
        multiple resources associated with it.

        Returns
        -------
        int

        See Also
        --------
        __len__
        """
        return sum(1 for _ in self.identifiers())

    def __len__(self) -> int:
        """Get the length of the object, which corresponds to the number of namespaces.

        Returns
        -------
        int

        See Also
        --------
        number_of_resources
        """
        return sum(1 for _ in self)

    def delete_annotation(self, value: str) -> None:
        """Delete a resource based on its identifier value.

        Parameters
        ----------
        value: str
            Identifier value of the resource to delete.

        Raises
        ------
        ValueError
            If no resource was found for the provided value.
        """
        for _, entries in self.objects():
            for entry in entries:
                if entry.identifier == value:
                    entry.remove_from_parent()
                    return
        raise ValueError(f"No annotation found for '{value}'")

    def to_dict(self) -> Dict[str, List[str]]:
        """Convert the simplified annotations to a dictionary.

        The keys of the dictionary represent namespaces and the values are lists of
        corresponding identifiers.

        Returns
        -------
        dict
        """
        return dict(self)

    def copy(self) -> Dict[str, List[str]]:
        """Convert the simplified annotations to a dictionary.

        This method returns a dictionary and not a new `SimplifiedAnnotationInterface`,
        since the interface is in place for compatibility purposes. Existing code will
        expect a copy of an object's `annotation` attribute to be a dictionary. This
        method is thus equivalent to `to_dict`.

        Returns
        -------
        dict

        See Also
        --------
        to_dict
        """
        # Example of usage of this method:
        # https://github.com/draeger-lab/MassChargeCuration/blob/main/MCC/ModelInterface/CobraPyInterface.py
        return self.to_dict()

    def clear(self) -> None:
        """Remove all annotations."""
        # Could probably handle this better by directly removing all standardized
        # annotations. Currently, empty StandardizedAnnotation objects will remain.
        for k in self:
            del self[k]

    def __deepcopy__(self, memo: dict):
        """Create a deepcopy of the SimplifiedAnnotationStore.

        Parameters
        ----------
        memo: dict
            Automatically passed parameter, dict of already copied items.

        Returns
        -------
        SimplifiedAnnotationInterface
            A new interface that is a deep copy of the original.
        """
        new = SimplifiedAnnotationInterface(self._metadata.__deepcopy__(memo))
        # memo[id(self)] = new
        return new
