"""Classes and functions to handle standardized annotation resources and qualifiers.

Resources represent a piece of information that relates to a cobrapy object through a
Qualifier, as captured in the  StandardizedAnnotation class.
"""

import logging
import re
from enum import Enum
from typing import (
    Any,
    Dict,
    Optional,
    Tuple,
    Union,
)


LOGGER = logging.getLogger(__name__)

__all__ = ["URL_IDENTIFIERS_PATTERN", "parse_identifiers_uri"]

# the URL pattern to parse namespace and identifier
URL_IDENTIFIERS_PATTERN = re.compile(r"^https?://identifiers.org/(.+?)[:/](.+)")


class Qualifier(Enum):
    """The possible qualifiers of a standardized annotation.

    The qualifiers and their detailed description are present in:
    https://co.mbine.org/author/biomodels.net-qualifiers/

    Qualifiers are divided into two groups

    Biological (bqb/bqbiol): These kinds of qualifiers define the relationship between a
    biological object represented by a model element and its annotation.

    Modelling (bqm/bqmodel): These kinds of qualifiers define the relationship between a
    modelling object and its annotation.

    See Also
    --------
    StandardizedAnnotation
    Resource
    """

    Biological_is = "bqb_is"
    """The biological entity represented by the cobrapy object is the subject of the
    referenced resource. This could serve to link a reaction to its counterpart in
    (e.g.) the ChEBI or Reactome databases."""
    Biological_hasPart = "bqb_hasPart"
    """The biological entity represented by the cobrapy object includes the subject of
    the referenced resource, either physically or logically. This relation might be used
    to link a complex to a description of its components"""
    Biological_isPartOf = "bqb_isPartOf"
    """The biological entity represented by the cobrapy object is a physical or logical
    part of the subject of the referenced resource. This relation might be used to link
    a component to the description of the complex to which it belongs."""
    Biological_isVersionOf = "bqb_isVersionOf"
    """The biological entity represented by the cobrapy object is a version or an
    instance of the subject of the referenced resource. This relation can be used to
    link a reaction to its EC-code."""
    Biological_hasVersion = "bqb_hasVersion"
    """The subject of the referenced resource is a version or an instance of the
    biological entity represented by the cobrapy object. This relation may be used to
    describe the components of a pool of metabolites."""
    Biological_isHomologTo = "bqb_isHomologTo"
    """The biological entity represented by the cobrapy object is homolog to the subject
    of the referenced resource, i.e. they share a common ancestor."""
    Biological_isDescribedBy = "bqb_isDescribedBy"
    """The biological entity represented by the cobrapy object is described by the
    referenced resource. This relation could be used, for example, to link a species or
    a parameter to a publication describing the quantity of the species or the value of
    the parameter."""
    Biological_isEncodedBy = "bqb_isEncodedBy"
    """The biological entity represented by the cobrapy object is encoded, either
    directly or by virtue of transitivity, by the subject of the referenced resource."""
    Biological_encodes = "bqb_encodes"
    """The biological entity represented by the cobrapy object encodes, either directly
    or by virtue of transitivity, the subject of the referenced resource."""
    Biological_occursIn = "bqb_occursIn"
    """The biological entity represented by the cobrapy object takes place in the
    subject of the reference resource"""
    Biological_hasProperty = "bqb_hasProperty"
    """The subject of the referenced resource is a property of the biological entity
    represented by the cobrapy object. This relation might be used when a biological
    entity has a given activity or exerts a specific function."""
    Biological_isPropertyOf = "bqb_isPropertyOf"
    """The biological entity represented by the cobrapy object is a property of the
    referenced resource."""
    Biological_hasTaxon = "bqb_hasTaxon"
    """The biological entity represented by the cobrapy object is taxonomically
    restricted, where the restriction is the subject of the referenced resource. This
    relation may be used to ascribe a species restriction to a biochemical reaction."""
    Biological_unknown = "bqb_unknown"
    """The relation between the resource an biological entity represented by the cobrapy
    object is unknown. Use sparingly, since this qualifier does not provide much
    information."""

    Modelling_is = "bqm_is"
    """The cobrapy object is the subject of the referenced resource. This may, for
    example, be used to link the model to an entry in a model database."""
    Modelling_isDescribedBy = "bqm_isDescribedBy"
    """The cobrapy object is described by the referenced resource. This could link a
    component (e.g., a reaction) to a publication describing it."""
    Modelling_isDerivedFrom = "bqm_isDerivedFrom"
    """The cobrapy object is derived from the modeling object represented by the
    referenced resource. For instance, they can be the fruit of a refinement or their
    adaptation for use in a different context."""
    Modelling_isInstanceOf = "bqm_isInstanceOf"
    """The cobrapy object is an instance of the subject of the referenced resource."""
    Modelling_hasInstance = "bqm_hasInstance"
    """The subject of the referenced resource is an instance of the cobrapy object. This
    could be used, for example, to link a generic model with its specific forms."""
    Modelling_unknown = "bqm_unknown"
    """The relation between the cobrapy object an the resource is unknown. Use
    sparingly, since this qualifier does not provide much information."""


class QualifiersAlias(set, Enum):
    Any_is = {Qualifier.Biological_is, Qualifier.Modelling_is}
    Any_isDescribedBy = {
        Qualifier.Biological_isDescribedBy,
        Qualifier.Modelling_isDescribedBy,
    }
    Any_unknown = {Qualifier.Biological_unknown, Qualifier.Modelling_unknown}
    Roughly_equals = {
        Qualifier.Biological_is,
        Qualifier.Biological_encodes,
        Qualifier.Biological_isEncodedBy,
        Qualifier.Biological_isVersionOf,
        Qualifier.Modelling_isInstanceOf,
    }
    Biological_any = {
        Qualifier.Biological_is,
        Qualifier.Biological_hasPart,
        Qualifier.Biological_isPartOf,
        Qualifier.Biological_isVersionOf,
        Qualifier.Biological_hasVersion,
        Qualifier.Biological_isHomologTo,
        Qualifier.Biological_isEncodedBy,
        Qualifier.Biological_encodes,
        Qualifier.Biological_isDescribedBy,
        Qualifier.Biological_hasTaxon,
        Qualifier.Biological_hasProperty,
        Qualifier.Biological_isPropertyOf,
        Qualifier.Biological_occursIn,
        Qualifier.Biological_unknown,
    }
    Biological_known = {
        Qualifier.Biological_is,
        Qualifier.Biological_hasPart,
        Qualifier.Biological_isPartOf,
        Qualifier.Biological_isVersionOf,
        Qualifier.Biological_hasVersion,
        Qualifier.Biological_isHomologTo,
        Qualifier.Biological_isEncodedBy,
        Qualifier.Biological_encodes,
        Qualifier.Biological_isDescribedBy,
        Qualifier.Biological_hasTaxon,
        Qualifier.Biological_hasProperty,
        Qualifier.Biological_isPropertyOf,
        Qualifier.Biological_occursIn,
    }
    Modelling_any = {
        Qualifier.Modelling_is,
        Qualifier.Modelling_isInstanceOf,
        Qualifier.Modelling_hasInstance,
        Qualifier.Modelling_isDescribedBy,
        Qualifier.Modelling_isDerivedFrom,
        Qualifier.Modelling_unknown,
    }
    Modelling_known = {
        Qualifier.Modelling_is,
        Qualifier.Modelling_isInstanceOf,
        Qualifier.Modelling_hasInstance,
        Qualifier.Modelling_isDescribedBy,
        Qualifier.Modelling_isDerivedFrom,
    }
    Roughly_instanceOf = {
        Qualifier.Modelling_isInstanceOf,
        Qualifier.Modelling_isDerivedFrom,
        Qualifier.Biological_isVersionOf,
    }


class Resource:
    """Defines a standardized annotation resource.

    Together with a Qualifier, Resource objects form the basis of a
    StandardizedAnnotation object. Resources are based around perennial URIs that link
    to scientific information and identifiers. A URI provided to the Resource class
    should be an https://identifiers.org URL.
    """

    def __init__(self, uri: str, strict: bool = True) -> None:
        """Initialize a standardized annotation Resource from a URI.

        Parameters
        ----------
        uri: str
            URI to use to create the resource. It should be of the format
            http(s)://identifiers.org/<namespace>/<identifier>.
        strict: bool, optional
            Whether to raise a ValueError when the provided URI does not match the
            identifiers.org pattern. If set to False, it will accept any URI and set
            `namespace` and `identifier` to None if the URI cannot be parsed.
            Default True.

        Raises
        ------
        ValueError
            If `strict` is set to True and a provided URI does not match the
            identifiers.org pattern.
        """
        self._namespace = None
        self._identifier = None
        self._parent = None
        self._strict = strict

        self.uri = uri

    def _set_parent(self, parent):
        if self._parent is None or self._parent is parent:
            self._parent = parent
        else:
            raise ValueError(
                "Resource already has a different parent. Create a new "
                "resource if you would like to add a resource to a second object."
            )

    def remove_from_parent(self):
        """Remove the Resource from its parent StandardizedAnnotation object."""
        if self._parent is None:
            raise ValueError(
                "Cannot remove resource, since no parent is associated with resource."
            )
        self._parent._remove_resource(self)
        self._parent = None

    @classmethod
    def from_data(
        cls,
        data: Union[Dict[str, str], Tuple[str, str], str, "Resource"],
        strict: bool = True,
    ):
        """Create a Resource instance from various formats of data.

        Parameters
        ----------
        data: dict, tuple, str or Resource
            If `data` is a dictionary and has a 'uri' key, the corresponding value will
            be used to create a Resource object, otherwise the keys 'namespace' and
            'identifier' will be converted an identifiers.org URI. If `data` is a tuple,
            it should be of length 2, where the first element is the namespace and the
            second the identifier. If `data` is a string, it is intrepreted as the URI
            of the Resource. If a Resource object is provided, this object will simply
            be returned.
        strict: bool, optional
            Whether to raise a ValueError when the provided URI does not match the
            identifiers.org pattern. If set to False, it will accept any URI and set
            `namespace` and `identifier` to None if the URI cannot be parsed.
            Default True.

        Returns
        -------
        Resource

        Raises
        ------
        TypeError
            If `data` is not of the correct type to creata a Resource object.
        ValueError
            If `strict` is set to True and a provided URI does not match the
            identifiers.org pattern.
        """
        if isinstance(data, Resource):
            data._strict = strict
            return data
        if isinstance(data, str):
            return Resource(data, strict=strict)
        if not isinstance(data, (dict, tuple)):
            raise TypeError(
                f"Resources can be created from str, tuple, dict, or Resource"
                f"types, not {type(data)}."
            )
        if isinstance(data, dict):
            if (uri := data.get("uri", None)) is not None:
                return Resource(uri, strict=strict)
            namespace = data["namespace"].lower()
            identifier = data["identifier"]
        else:
            namespace, identifier = data
        if not isinstance(namespace, str) or not isinstance(identifier, str):
            raise TypeError("Namespace and identifier should be of type str.")
        uri = f"https://identifiers.org/{namespace}/{identifier}"
        return Resource(uri, strict=strict)

    @property
    def uri(self) -> str:
        """Get the URI of the resource.

        Returns
        -------
        str
        """
        return self._uri

    @uri.setter
    def uri(self, value: str) -> None:
        """Set the URI of the resource.

        Parameters
        ----------
        value: str
            The URI, typically of the format
            'https://identifiers.org/<namespace>/<identifier>'.
        """
        if re.match(URL_IDENTIFIERS_PATTERN, value):
            identifier_match = parse_identifiers_uri(value)

            if identifier_match is not None:
                namespace, identifier = identifier_match
                self._namespace = namespace
                self._identifier = identifier
                self._uri = value
                return
        if self._strict:
            raise ValueError(
                f"The provided URI is not a valid identifiers.org address: {value}"
            )
        self._namespace = None
        self._identifier = None
        self._uri = value

    @property
    def namespace(self) -> Optional[str]:
        """Get the namespace of the resource.

        Returns
        -------
        str or None
            The namespace of the resource (e.g. 'chebi'). If `strict` was set to False
            upon initialization and the URI could not be interpreted, this method will
            return None.
        """
        return self._namespace

    @property
    def identifier(self) -> Optional[str]:
        """Get the identifier of the resource.

        Returns
        -------
        str or None
            The identifier of the resource (e.g. 'CHEBI:36927'). If `strict` was set to
            False upon initialization and the URI could not be interpreted, this method
            will return None.
        """

        return self._identifier

    def to_dict(self):
        """Convert the resource to a dictionary.

        The dictionary will have the key "uri" and optionally "namespace" and
        "identifier", if the corresponding attributes are not None.

        Returns
        -------
        dict
        """
        return {
            k: v
            for k in ["uri", "namespace", "identifier"]
            if (v := getattr(self, k, None)) is not None
        }

    def __repr__(self) -> str:
        """Get the string representation of the resource.

        Returns
        -------
        str
        """
        return (
            f"{self.__class__.__module__}.{self.__class__.__qualname__}"
            f"({self.to_dict()})"
        )

    def _repr_html(self) -> str:
        """Return the resource as an HTML string.

        Returns
        -------
        str
        """
        return f"""
        <table>
            <tr>
                <td><strong>URI</strong></td><td>{self.uri}</td>
            </tr><tr>
                <td><strong>Namespace</strong></td><td>{self.namespace}</td>
            </tr><tr>
                <td><strong>Identifier</strong></td><td>{self.identifier}</td>
            </tr><tr>
                <td><strong>Memory address</strong></td>
                <td>{id(self):#x}</td>
            </tr>
        </table>"""

    def __eq__(self, other: Any):
        """Determine equality between the resource and another object.

        Parameters
        ----------
        other: Resource, str, dict or tuple
            If `other` is a Resource object, the objects are equal if their URIs are
            equal. Similarly, if `other` is a string or a dict with a 'uri' key, these
            values will be compared to the `uri` attribute of this object. If other is a
            tuple, the first element wil be compared to the `namespace` property of this
            object and the second element to the `identifier` property.

        Returns
        -------
        bool
            True if objects are equal, based on their URI (or namespace + identifier
            combination, in the case of comparison to a tuple). Returns False otherwise,
            including when `other` could not be interpreted as a Resource.
        """
        # TODO: Should we consider https://identifiers.org/... and
        # http://identifiers.org/... URIs as equal? Maybe we should convert everything
        # to https upon initialization.
        if isinstance(other, str):
            return self.uri == other
        if isinstance(other, Resource):
            return self.uri == other.uri
        if isinstance(other, dict):
            return self.uri == other.get("uri")
        if isinstance(other, tuple):
            if len(other) != 2 or self.namespace is None:
                return False
            return self.namespace == other[0] and self.identifier == other[1]
        return False

    def __hash__(self):
        """Create a hash based on the URI of the resource."""
        return hash(self.uri)


def parse_identifiers_uri(uri: str) -> Optional[Tuple[str, str]]:
    """Parse namespace and term from given identifiers annotation uri.

    Parameters
    ----------
    uri : str
        uri (identifiers.org url)

    Returns
    -------
    (namespace, identifier) if resolvable, None otherwise
    """
    match = URL_IDENTIFIERS_PATTERN.match(uri)
    if match:
        namespace, identifier = match.group(1), match.group(2)
        if namespace.isupper():
            identifier = f"{namespace}:{identifier}"
            namespace = namespace.lower()
    else:
        LOGGER.warning(
            f"{uri} does not conform to "
            f"'http(s)://identifiers.org/collection/id' or"
            f"'http(s)://identifiers.org/COLLECTION:id"
        )
        return None
    return namespace, identifier


DEFAULT_QUALIFIERS = {
    "pubmed": Qualifier.Modelling_isDescribedBy,
    "doi": Qualifier.Modelling_isDescribedBy,
    "ec-code": Qualifier.Biological_isVersionOf,
    "go": Qualifier.Biological_isVersionOf,
    "eco": Qualifier.Modelling_isDescribedBy,
    "google.patent": Qualifier.Modelling_isDescribedBy,
    "taxonomy": Qualifier.Biological_hasTaxon,
    "arxiv": Qualifier.Modelling_isDescribedBy,
    "isbn": Qualifier.Modelling_isDescribedBy,
    "bigg.model": Qualifier.Modelling_is,
}


def get_default_qualifier(namespace):
    return DEFAULT_QUALIFIERS.get(str(namespace).lower(), Qualifier.Biological_is)
    # return Qualifier.Biological_is
