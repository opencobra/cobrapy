"""
Provide a concrete implementation of the carveme repository interface.
"""


from io import BytesIO

import httpx

from .abstract_model_repository import AbstractModelRepository


def _decode_model_path(model_path):
    """Decode the model path to EMBL GEMs."""
    tokens = model_path.split("_")
    genus  = tokens[0]

    directory = genus.lower()
    alphabet  = directory[0]

    return f"{alphabet}/{directory}/{model_path}"

class EMBLGems(AbstractModelRepository):
    """
    Define a concrete implementation of the EMBL GEMs repository.

    Attributes
    ----------
    name : str
        The name of the EMBL GEMs repository.

    """

    name: str = "EMBL GEMs"

    def __init__(
        self,
        **kwargs,
    ) -> None:
        """
        Initialize a EMBL GEMs repository interface.

        Other Parameters
        ----------------
        kwargs
            Passed to the parent constructor in order to enable multiple inheritance.

        """
        super().__init__(url="https://github.com/cdanielmachado/embl_gems/blob/master/models/", **kwargs)

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

        decoded_path = _decode_model_path(model_id)

        filename = f"{model_id}.xml.gz"
        print(self._url.join(decoded_path).join(filename))
        with self._progress, httpx.stream(
            method="GET", url=self._url.join(decoded_path).join(filename),
            params={"raw": "true"},
            follow_redirects=True
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
