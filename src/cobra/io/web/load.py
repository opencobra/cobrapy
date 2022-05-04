"""Provide a function ``load_model`` to access remote model repositories."""


import gzip
import logging
from typing import TYPE_CHECKING, Iterable

import diskcache
import httpx
import libsbml

from ...core import Configuration
from ..sbml import _sbml_to_model
from .abstract_model_repository import AbstractModelRepository
from .bigg_models_repository import BiGGModels
from .biomodels_repository import BioModels
from .cobrapy_repository import Cobrapy


if TYPE_CHECKING:
    from cobra.core import Model


logger = logging.getLogger(__name__)
configuration = Configuration()


DEFAULT_REPOSITORIES = (
    Cobrapy(),
    BiGGModels(),
    BioModels(),
)


def load_model(
    model_id: str,
    repositories: Iterable[AbstractModelRepository] = DEFAULT_REPOSITORIES,
    cache: bool = True,
) -> "Model":
    """
    Download an SBML model from a remote repository.

    Downloaded SBML documents are by default stored in a cache on disk such that future
    access is much faster. By default, models can be loaded from the following
    repositories:

    * BiGG Models
    * BioModels

    You can use the ``AbstractModelRepository`` class as a parent to implement your own
    repository accessor which you pass to the ``load_model`` function. In case you
    implement a new interface, please consider submitting a pull request to COBRApy.

    Parameters
    ----------
    model_id : str
        The identifier of the desired metabolic model. This is typically repository
        specific.
    repositories : iterable, optional
        An iterable of repository accessor instances. The model_id is searched in order.
    cache : bool, optional
        Whether or not to use the local caching mechanism (default yes).

    Returns
    -------
    Model
        A model instance generated from the SBML document.

    Raises
    ------
    RuntimeError
        As with any internet connection, there are multiple errors that can occur.

    Examples
    --------
    # Most of the time calling `load_model` with an identifier should be enough.
    >>> print(load_model("e_coli_core"))
    e_coli_core
    >>> print(load_model("MODEL1510010000"))
    MODEL1510010000

    See Also
    --------
    BiGGModels
    BioModels

    """
    if cache:
        data = _cached_load(
            model_id=model_id,
            repositories=repositories,
        )
    else:
        data = _fetch_model(model_id=model_id, repositories=repositories)
    return get_model_from_gzip_sbml(data)


def _cached_load(
    model_id: str,
    repositories: Iterable[AbstractModelRepository],
) -> bytes:
    """
    Attempt to load a gzip-compressed SBML document from the cache.

    If the given model identifier is not in the cache, the remote repositories are
    searched.

    Parameters
    ----------
    model_id : str
        The identifier of the desired metabolic model. This is typically repository
        specific.
    repositories : iterable
        An iterable of repository accessor instances. The model_id is searched in order.

    Returns
    -------
    bytes
        A gzip-compressed, UTF-8 encoded SBML document.

    """
    with diskcache.Cache(
        directory=str(configuration.cache_directory),
        size_limit=configuration.max_cache_size,
    ) as cache:
        try:
            return cache[model_id]
        except KeyError:
            data = _fetch_model(model_id=model_id, repositories=repositories)
            cache.set(key=model_id, value=data, expire=configuration.cache_expiration)
            return data


def _fetch_model(
    model_id: str,
    repositories: Iterable[AbstractModelRepository],
) -> bytes:
    """
    Attempt to load a gzip-compressed SBML document from the given repositories.

    Parameters
    ----------
    model_id : str
        The identifier of the desired metabolic model. This is typically repository
        specific.
    repositories : iterable
        An iterable of repository accessor instances. The model_id is searched in order.

    Returns
    -------
    bytes
        A gzip-compressed, UTF-8 encoded SBML document.

    """
    for repository in repositories:
        logger.info(
            f"Attempting to fetch '{model_id}' from the {repository.name} repository."
        )
        try:
            return repository.get_sbml(model_id=model_id)
        except OSError:
            logger.debug(
                f"Model '{model_id} not found in the local "
                f"repository {repository.name}.'"
            )
        except httpx.HTTPStatusError as error:
            if error.response.status_code == 404:
                logger.debug(
                    f"Model '{model_id}' not found in the {repository.name} repository."
                )
                continue
            raise RuntimeError(
                f"The connection to the {repository.name} repository failed."
            ) from error
        except httpx.RequestError as error:
            raise RuntimeError(
                f"The connection to the {repository.name} repository failed."
            ) from error
    raise RuntimeError(
        f"The model '{model_id}' could not be found in any of the repositories."
    )


def get_model_from_gzip_sbml(stream: bytes) -> "Model":
    """
    Generate a model instance from a gzip-compressed, UTF-8 encoded SBML document.

    Parameters
    ----------
    stream : bytes
        A gzip-compressed, UTF-8 encoded SBML document.

    Returns
    -------
    Model
        A model instance generated from the SBML document.

    """
    return _sbml_to_model(
        libsbml.readSBMLFromString(gzip.decompress(stream).decode("utf-8"))
    )
