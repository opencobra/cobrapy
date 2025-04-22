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
    Biological (bqb)    These kinds of qualifiers define the relationship between a
                        biological object represented by a model element and its
                        annotation.
    Modelling (bqm)     These kinds of qualifiers define the relationship between a
                        modelling object and its annotation.

    See Also
    --------
    StandardizedAnnotation
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


class Resource:
    def __init__(self, uri: str) -> None:
        self._namespace = None
        self._identifier = None
        self._parent = None

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
        if self._parent is None:
            raise ValueError(
                "Cannot remove resource, since no parent is associated with resource."
            )
        self._parent._remove_resource(self)
        self._parent = None

    @classmethod
    def from_data(cls, data: Union[Dict[str, str], Tuple[str, str], str, "Resource"]):
        if isinstance(data, Resource):
            return data
        if isinstance(data, str):
            return Resource(data)
        if not isinstance(data, (dict, tuple)):
            raise TypeError(
                f"Resources can be created from str, tuple, dict, or Resource"
                f"types, not {type(data)}."
            )
        if isinstance(data, dict):
            if (uri := data.get("uri", None)) is not None:
                return Resource(uri)
            namespace = data["namespace"].lower()
            identifier = data["identifier"]
        else:
            namespace, identifier = data
        if not isinstance(namespace, str) or not isinstance(identifier, str):
            raise TypeError("Namespace and identifier should be of type str.")
        uri = f"https://identifiers.org/{namespace}/{identifier}"
        return Resource(uri)

    @property
    def uri(self) -> str:
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
            self._uri = value
        else:
            # TODO: Warn user
            self._namespace = None
            self._identifier = None
            self._uri = value

    @property
    def namespace(self) -> Optional[str]:
        return self._namespace

    @property
    def identifier(self) -> Optional[str]:
        return self._identifier

    def to_dict(self):
        return {
            k: v
            for k in ["uri", "namespace", "identifier"]
            if (v := getattr(self, k, None)) is not None
        }

    def __repr__(self) -> str:
        return (
            f"{self.__class__.__module__}.{self.__class__.__qualname__}"
            f"({self.to_dict()})"
        )

    def _repr_html(self) -> str:
        """Return the resource as an HTML string."""
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
