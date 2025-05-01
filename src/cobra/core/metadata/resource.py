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

__all__ = ["parse_identifiers_uri", "namespace_and_identifier_to_uri"]

# the URL pattern to parse namespace and identifier
URL_OLD_IDENTIFIERS_PATTERN = re.compile(
    r"^https?://identifiers.org/(.+?)/(.+)(\?.*)?$"
)
URL_COMPACT_IDENTIFIERS_PATTERN = re.compile(
    r"^https?://identifiers.org/((.+)/)?([^:/]+):(.+)(\?.*)?$"
)
COMPACT_IDENTIFIERS_PATTERN = re.compile(r"^[a-zA-Z0-9_\.\-]+:.+$")
# Code to update the following constants can be found in:
# scripts/parse_identifiers_registry.py
COMPACT_URL_INTEGRATED_NAMESPACES = [
    "mge",
    "ark",
    "bto",
    "cco",
    "cl",
    "chebi",
    "cheminf",
    "did",
    "envo",
    "eco",
    "fma",
    "foodon",
    "gsso",
    "go",
    "go_ref",
    "gro",
    "doid",
    "hp",
    "mir",
    "mp",
    "ms",
    "mcro",
    "mi",
    "ma",
    "mgi",
    "nando",
    "nmr",
    "oma.hog",
    "ocid",
    "opl",
    "obcs",
    "pw",
    "pato",
    "eo",
    "po",
    "mod",
    "pr",
    "rrid",
    "so",
    "swh",
    "stato",
    "slm",
    "sbo",
    "uberon",
    "uo",
    "mzspec",
    "vario",
]
COMPACT_URL_IDENTIFIERS_WITH_COLON = [
    "arraymap",
    "bgee.family",
    "bgee.organ",
    "bgee.stage",
    "biocyc",
    "bbkg",
    "cabri",
    "dip",
    "ga4ghdos",
    "dev.ga4ghdos",
    "doi",
    "glyconavi",
    "gramene.gene",
    "gramene.taxonomy",
    "hgnc",
    "imgt.hla",
    "isbn",
    "kegg.environ",
    "kegg.genes",
    "kegg",
    "metacyc.compound",
    "metacyc.reaction",
    "miriam.collection",
    "miriam.resource",
    "narcis",
    "nbn",
    "nmdc",
    "orphanet.ordo",
    "panther.family",
    "ps",
    "psipar",
    "sisu",
    "storedb",
    "tair.gene",
    "tair.protein",
    "treebase",
    "vgnc",
]
# This does not actually fix anything, since these special cases are broken on
# identifiers.org. However, future cases could need this logic.
COMPACT_URL_NAMESPACE_EXCEPTIONS = {
    "hog": "oma.hog",  # This one is weird. Hard to check what works, server is down.
    "peo": "eo",  # The sampel URL on identifiers.org is wrong and does not work.
}


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
    QualifiersAlias
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
    """Aliases for common combinations of qualifiers.

    Examples
    --------
    >>> from cobra.core.metadata import Qualifier, QualifiersAlias
    >>> Qualifier.Biological_isPartOf in QualifiersAlias.Biological_any
    True
    >>> Qualifier.Modelling_unknown in QualifiersAlias.Modelling_known
    False

    See Also
    --------
    Qualifier
    StandardizedAnnotationStore.resources_for
    """

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

    Parameters
    ----------
    uri: str
        URI to use to create the resource. It should be of the old identifiers.org
        format 'http(s)://identifiers.org/<namespace>/<identifier>' or the compact
        identifier URL format
        http(s)://identifiers.org/<namespace>:<identifier>. Alternatively, a
        compact identifier can be provided directly, without the preceding
        'http(s)://identifiers.org/', e.g. 'CHEBI:11881'.
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

    Examples
    --------
    >>> from cobra.core import Resource
    >>> chebi_resource = Resource("https://identifiers.org/chebi/CHEBI:11881")
    >>> chebi_resource.namespace
    chebi
    >>> chebi_resource.identifier
    CHEBI:11881
    >>> chebi_resource == Resource("CHEBI:11881")
    True
    >>> try:
            ebi_resource = Resource(
                "https://www.ebi.ac.uk/chebi/searchId.do?chebiId=CHEBI:11881",
            )
        except ValueError:
            print("Not an identifiers.org URL.")
    Not an identifiers.org URL.
    >>> ebi_resource = Resource(
            "https://www.ebi.ac.uk/chebi/searchId.do?chebiId=CHEBI:11881",
            strict=False,
        )
    >>> ebi_resource.uri
    https://www.ebi.ac.uk/chebi/searchId.do?chebiId=CHEBI:11881
    >>> ebi_resource.namespace is None
    True
    >>> chebi_resource == ebi_resource
    False
    """

    def __init__(self, uri: str, strict: bool = True) -> None:
        """Initialize a standardized annotation Resource from a URI."""
        self._namespace = None
        self._identifier = None
        self._parent = None
        self._strict = strict

        self.uri = uri

    def _set_parent(self, parent):
        if self._parent is None or parent is None:
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
            second the identifier. If `data` is a string, it is interpreted as the URI
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
            If `data` is not of the correct type to create a Resource object.
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
            if len(data) != 2:
                raise TypeError(
                    "Tuples should have length 2 to be converted to a Resource, "
                    f"not length {len(data)}: {data}"
                )
            namespace = data[0].lower()
            identifier = data[1]
        if not isinstance(namespace, str) or not isinstance(identifier, str):
            raise TypeError("Namespace and identifier should be of type str.")
        uri = namespace_and_identifier_to_uri(namespace, identifier)
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
            The URI or identifiers.org compact identifier, typically of the format
            'https://identifiers.org/<namespace>:<identifier>'.
        """
        if (identifier_match := parse_identifiers_uri(value)) is not None:
            namespace, identifier, _provider, uri = identifier_match
            self._namespace = namespace
            self._identifier = identifier
            self._uri = uri
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
        return f"Resource({self.uri})"

    def _repr_html_(self) -> str:
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
            The value of `other` is converted to a Resource object, if it was not
            already of that type. If the Resource object has a `namespace` and
            `identifier` that are not None, i.e. the URI was a valid identifiers.org
            URI, the two objects are equal if their `namespace` and `identifier` are
            equal. If the `namespace` and `identifier` are None, the resources are equal
            if their URI is equal.
        Returns
        -------
        bool
            True if objects are equal, based on their `namespace` and `identifier`
            combination (or URI, if those attributes are None).
        """
        if isinstance(other, (dict, tuple, str)):
            try:
                other_resource = Resource.from_data(other, strict=False)
            except TypeError:
                return False
            return self == other_resource
        if isinstance(other, Resource):
            if (
                self.namespace is None
                or self.identifier is None
                or other.namespace is None
                or other.identifier is None
            ):
                return self.uri == other.uri
            return (
                self.namespace == other.namespace
                and self.identifier == other.identifier
            )
        return False

    def __hash__(self):
        """Create a hash of the namespace and identifier, or URI of the resource."""
        if self.namespace is None:
            return hash(self.uri)
        else:
            return hash((self.namespace, self.identifier))


def parse_identifiers_uri(uri: str) -> Optional[Tuple[str, str, Optional[str], str]]:
    """Parse namespace and term from given identifiers annotation uri.

    Parameters
    ----------
    uri : str
        uri (identifiers.org url) or identifiers.org compact identifier (e.g.
        "CHEBI:11881").

    Returns
    -------
    (namespace, identifier, provider) if resolvable, None otherwise
    """
    if not (uri.startswith("http://") or uri.startswith("https://")):
        # Try to interpret the uri as a identifiers.org compact identifier.
        if not COMPACT_IDENTIFIERS_PATTERN.match(uri):
            return None
        uri = f"https://identifiers.org/{uri}"
    # Try to match the new format first
    match = URL_COMPACT_IDENTIFIERS_PATTERN.match(uri)
    if match:
        provider, orig_namespace, identifier = (
            match.group(2),
            match.group(3),
            match.group(4),
        )
        # For most compact URLs the namespace prefix is simply the lower case
        # version of what is matched in the URL, but there are some exceptions
        # that need correcting.
        namespace = COMPACT_URL_NAMESPACE_EXCEPTIONS.get(
            orig_namespace.lower(), orig_namespace.lower()
        )
        # If what is interpreted as provider is a namespace where identifiers have
        # colons, the url should be intrepreted as the old format instead.
        if provider == namespace or provider not in COMPACT_URL_IDENTIFIERS_WITH_COLON:
            # In the cases where the namespace is integrated in the compact URL, the
            # identifier should be reconstructed. E.g. in the case of ChEBI, the
            # namespace is 'chebi' and a identifier can be 'CHEBI:11881'.
            if namespace in COMPACT_URL_INTEGRATED_NAMESPACES:
                identifier = f"{orig_namespace}:{identifier}"
            return namespace, identifier, provider, uri
    # Otherwise try the old format
    match = URL_OLD_IDENTIFIERS_PATTERN.match(uri)
    if match:
        provider, namespace, identifier = (
            None,
            match.group(1),
            match.group(2),
        )
        return namespace, identifier, provider, uri

    LOGGER.warning(
        f"{uri} does not conform to "
        f"'http(s)://identifiers.org/namespace/id' or "
        f"'http(s)://identifiers.org/(provider/)namespace:id"
    )
    return None


def namespace_and_identifier_to_uri(namespace: str, identifier: str) -> str:
    """Convert a namespace and identifier pair to a identifiers.org URL.

    Parameters
    ----------
    namespace : str
        Namespace of the entity.
    identifier : str
        Identifier of the entity.

    Returns
    -------
    str
        identifiers.org URL
    """
    if namespace in COMPACT_URL_INTEGRATED_NAMESPACES:
        return f"https://identifiers.org/{namespace}/{identifier}"
    return f"https://identifiers.org/{namespace}:{identifier}"


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
    return DEFAULT_QUALIFIERS.get(str(namespace), Qualifier.Biological_is)
    # return Qualifier.Biological_is
