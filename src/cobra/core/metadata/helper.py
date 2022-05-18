"""Helper functions and properties for dealing with annotation data.

This is a unification of helper functions from sbml.py and cvterm.py.

"""
import logging
import re

from typing import Union, Tuple

LOGGER = logging.getLogger(__name__)

__all__ = ['URL_IDENTIFIERS_PATTERN', '_parse_identifiers_uri']

# the URL pattern to parse provider and identifier
URL_IDENTIFIERS_PATTERN = re.compile(r"^https?://identifiers.org/(.+?)[:/](.+)")


def _parse_identifiers_uri(uri: str) -> Union[None, Tuple[str, str]]:
    """Parse provider and term from given identifiers annotation uri.

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
            identifier = f"{provider}:{identifier}"
            provider = provider.lower()
    else:
        LOGGER.warning(
            f"{uri} does not conform to "
            f"'http(s)://identifiers.org/collection/id' or"
            f"'http(s)://identifiers.org/COLLECTION:id"
        )
        return None

    return provider, identifier
