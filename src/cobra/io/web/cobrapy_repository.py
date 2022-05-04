"""Provide functions for loading metabolic models from local package data."""


import importlib_resources

import cobra.data

from .abstract_model_repository import AbstractModelRepository


class Cobrapy(AbstractModelRepository):
    """
    Define a concrete implementation of the cobrapy (local package) repository.

     Attributes
    ----------
    name : str
        The name of the Cobrapy repository.
    """

    name: str = "Cobrapy"

    def __init__(
        self,
        **kwargs,
    ) -> None:
        """
        Initialize a local Cobrapy repository interface.

        Other Parameters
        ----------------
        kwargs
            Passed to the parent constructor in order to enable multiple inheritance.

        """
        super().__init__(url="file:////", **kwargs)

    def get_sbml(self, model_id: str) -> bytes:
        """
        Attempt to open an SBML document from the local repository.

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
        return (
            importlib_resources.files(cobra.data)
            .joinpath(f"{model_id}.xml.gz")
            .read_bytes()
        )
