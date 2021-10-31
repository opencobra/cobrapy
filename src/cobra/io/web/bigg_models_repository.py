"""Provide a concrete implementation of the BioModels repository interface."""


from io import BytesIO

import httpx

from .abstract_model_repository import AbstractModelRepository


class BiGGModels(AbstractModelRepository):
    """
    Define a concrete implementation of the BiGG Models repository.

    Attributes
    ----------
    name : str
        The name of the BiGG Models repository.

    """

    name: str = "BiGG Models"

    def __init__(
        self,
        **kwargs,
    ) -> None:
        """
        Initialize a BiGG Models repository interface.

        Other Parameters
        ----------------
        kwargs
            Passed to the parent constructor in order to enable multiple inheritance.

        """
        super().__init__(url="http://bigg.ucsd.edu/static/models/", **kwargs)

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
        compressed = BytesIO()
        filename = f"{model_id}.xml.gz"
        with self._progress, httpx.stream(
            method="GET", url=self._url.join(filename)
        ) as response:
            response.raise_for_status()
            task_id = self._progress.add_task(
                description="download",
                total=int(response.headers["Content-Length"]),
                model_id=model_id,
            )
            for chunk in response.iter_bytes():
                compressed.write(chunk)
                self._progress.update(task_id=task_id, advance=len(chunk))
        compressed.seek(0)
        return compressed.read()
