"""Provide functions for loading metabolic models over the wire."""


import gzip
from io import BytesIO
from typing import List

import httpx
import pydantic

from .abstract_model_repository import AbstractModelRepository


class BioModelsFile(pydantic.BaseModel):
    """Define a single BioModels file description."""

    name: str
    size: int = pydantic.Field(alias="fileSize")


class BioModelsFilesResponse(pydantic.BaseModel):
    """Define the BioModels files JSON response."""

    main: List[BioModelsFile] = []


class BioModels(AbstractModelRepository):
    """
    Define a concrete implementation of the BioModels repository.

    Attributes
    ----------
    name : str
        The name of the BioModels repository.

    """

    name: str = "BioModels"

    def __init__(
        self,
        **kwargs,
    ) -> None:
        """
        Initialize a BioModels repository interface.

        Other Parameters
        ----------------
        kwargs
            Passed to the parent constructor in order to enable multiple inheritance.

        """
        super().__init__(url="https://www.ebi.ac.uk/biomodels/model/", **kwargs)

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

        Raises
        ------
        httpx.HTTPError
            In case there are any connection problems.

        """
        data = BytesIO()
        response = httpx.get(
            url=self._url.join(f"files/{model_id}"),
            headers={"Accept": "application/json"},
        )
        response.raise_for_status()
        files = BioModelsFilesResponse.parse_obj(response.json())
        for model in files.main:
            if model.name.endswith("xml"):
                break
        else:
            RuntimeError(f"Could not find an SBML document for '{model_id}'.")
        with self._progress, httpx.stream(
            method="GET",
            url=self._url.join(f"download/{model_id}"),
            params={"filename": model.name},
        ) as response:
            response.raise_for_status()
            task_id = self._progress.add_task(
                description="download",
                total=model.size,
                model_id=model_id,
            )
            for chunk in response.iter_bytes():
                data.write(chunk)
                self._progress.update(task_id=task_id, advance=len(chunk))
        data.seek(0)
        return gzip.compress(data.read())
