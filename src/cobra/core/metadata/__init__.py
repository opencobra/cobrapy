"""Provide functions for loading, saving and modifying metadata annotations."""

from cobra.core.metadata.cvterm import CVTerm, CVTermList, ExternalResources, Qualifier
from cobra.core.metadata.helper import URL_IDENTIFIERS_PATTERN, parse_identifiers_uri
from cobra.core.metadata.history import Creator, History, HistoryDatetime
from cobra.core.metadata.keyvaluepairs import KeyValuePairs
from cobra.core.metadata.metadata import MetaData
