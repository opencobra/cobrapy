"""Provide an abstract base class that describes a remote model repository."""


from abc import ABC, abstractmethod
from typing import Union

import httpx
from rich.progress import (
    BarColumn,
    DownloadColumn,
    Progress,
    TextColumn,
    TimeRemainingColumn,
    TransferSpeedColumn,
)


class AbstractModelRepository(ABC):
    """
    Define an abstract base class that describes a remote model repository.

    Attributes
    ----------
    name : str
        The name of the remote repository.

    """

    _progress = Progress(
        TextColumn("{task.fields[model_id]}", justify="right"),
        BarColumn(bar_width=None),
        "[progress.percentage]{task.percentage:>3.1f}%",
        DownloadColumn(),
        TransferSpeedColumn(),
        TimeRemainingColumn(),
    )
    name: str = "Abstract"

    def __init__(self, *, url: Union[httpx.URL, str], **kwargs) -> None:
        """
        Initialize the model repository.

        Parameters
        ----------
        url : httpx.URL or str
            The base URL from where to load the models.

        Other Parameters
        ----------------
        kwargs
            Passed to the parent constructor in order to enable multiple inheritance.

        """
        super().__init__(**kwargs)
        self._url = httpx.URL(url=url)

    @property
    def url(self) -> httpx.URL:
        """Return the repository's URL."""
        return self._url.copy_with()

    @abstractmethod
    def get_sbml(self, model_id: str) -> bytes:
        """
        Attempt to download an SBML document from the repository.

        Parameters
        ----------
        model_id : str
            The identifier of the desired metabolic model. This is typically repository
            specific.

        Returns
        -------
        bytes
            A gzip-compressed, UTF-8 encoded SBML document.

        """
        raise NotImplementedError("Implement `get_sbml` in a concrete child class.")
