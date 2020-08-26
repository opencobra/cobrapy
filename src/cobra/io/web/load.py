"""Provide a function ``load_model`` to access remote model repositories."""


import gzip
import logging
import pathlib
from typing import TYPE_CHECKING, Iterable, Optional, Union

import appdirs
import diskcache
import httpx
import libsbml

from ..sbml import _sbml_to_model
from .abstract_model_repository import AbstractModelRepository
from .bigg_models_repository import BiGGModels
from .biomodels_repository import BioModels


if TYPE_CHECKING:
    from cobra.core import Model


logger = logging.getLogger(__name__)


def load_model(
    model_id: str,
    repositories: Iterable[AbstractModelRepository] = (BiGGModels(), BioModels()),
    cache: bool = True,
    cache_directory: Optional[Union[pathlib.Path, str]] = None,
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
    cache_directory : pathlib.Path or str, optional
        A path where the cache should reside if caching is desired. The default
        directory depends on the operating system.

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
        if cache_directory is None:
            cache_directory = pathlib.Path(
                appdirs.user_cache_dir(appname="cobrapy", appauthor="opencobra")
            )
        data = _cached_load(
            model_id=model_id,
            repositories=repositories,
            cache_directory=pathlib.Path(cache_directory),
        )
    else:
        data = _fetch_model(model_id=model_id, repositories=repositories)
    return get_model_from_gzip_sbml(data)


def _cached_load(
    model_id: str,
    repositories: Iterable[AbstractModelRepository],
    cache_directory: pathlib.Path,
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
    cache_directory : pathlib.Path
        A path where the cache should reside.

    Returns
    -------
    bytes
        A gzip-compressed, UTF-8 encoded SBML document.

    """
    if not cache_directory.is_dir():
        logger.debug(f"Creating cache directory '{str(cache_directory)}'.")
        cache_directory.mkdir(parents=True)
    with diskcache.Cache(directory=str(cache_directory)) as cache:
        try:
            return cache[model_id]
        except KeyError:
            data = _fetch_model(model_id=model_id, repositories=repositories)
            # (Midnighter): We could expire models after some time here.
            cache.set(key=model_id, value=data)
            return data


def _fetch_model(
    model_id: str, repositories: Iterable[AbstractModelRepository],
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
