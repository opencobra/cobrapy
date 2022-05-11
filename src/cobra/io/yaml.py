"""Provide functions for I/O in YAML format."""

import io
from pathlib import Path
from typing import TYPE_CHECKING, Any, Dict, Optional, Union

from ruamel.yaml.compat import StringIO
from ruamel.yaml.main import YAML

from .dict import model_from_dict, model_to_dict


if TYPE_CHECKING:
    from io import TextIOBase

    from cobra import Model


YAML_SPEC = "1.2"


class CobraYAML(YAML):
    """Define custom subclass for YAML I/O."""

    def dump(
        self, data: Dict, stream: Optional["TextIOBase"] = None, **kwargs: Any
    ) -> str:
        """Dump YAML data.

        Parameters
        ----------
        data : dict
            A dictionary representing the cobra model and its components.
        stream : TextIOBase, optional
            A text stream object inheriting from `io.TextIOBase`. If None,
            `ruamel.yaml.compat.StringIO` is used (default None).
        **kwargs : Any
            Keyword arguments passed on to `ruamel.yaml.main.YAML.dump`.

        Returns
        -------
        str
            YAML string representation of `data`.

        """
        if stream is None:
            inefficient = True
            stream = StringIO()
        else:
            inefficient = False
        YAML.dump(self, data, stream, **kwargs)
        if inefficient:
            return stream.getvalue()


yaml = CobraYAML(typ="rt")


def to_yaml(model: "Model", sort: bool = False, **kwargs: Any) -> str:
    """Return the model as a YAML string.

    Parameters
    ----------
    model : cobra.Model
        The cobra model to represent.
    sort : bool, optional
        Whether to sort the metabolites, reactions, and genes or maintain the
        order defined in the model (default False).
    **kwargs : Any
        Keyword arguments passed on to `CobraYAML.dump`.

    Returns
    -------
    str
        YAML string representation of the cobra model.

    See Also
    --------
    save_yaml_model : Write directly to a file.
    ruamel.yaml.dump : Base function.

    """
    obj = model_to_dict(model, sort=sort)
    obj["version"] = YAML_SPEC
    return yaml.dump(obj, **kwargs)


def from_yaml(document: str) -> "Model":
    """Load a cobra model from a YAML string.

    Parameters
    ----------
    document : str
        The YAML string representation of a cobra model.

    Returns
    -------
    cobra.Model
        The cobra model as interpreted from the YAML document.

    See Also
    --------
    load_yaml_model : Load directly from a file.

    """
    content = StringIO(document)
    return model_from_dict(yaml.load(content))


def save_yaml_model(
    model: "Model", filename: str, sort: bool = False, **kwargs: Any
) -> None:
    """Write the cobra model to a file in YAML format.

    Parameters
    ----------
    model : cobra.Model
        The cobra model to represent.
    filename : str or file-like
        File path or descriptor that the YAML representation should be
        written to.
    sort : bool, optional
        Whether to sort the metabolites, reactions, and genes or maintain the
        order defined in the model (default False).
    **kwargs : Any
        Keyword arguments passed to `CobraYAML.dump`.

    See Also
    --------
    to_yaml : Return a string representation.
    ruamel.yaml.dump : Base function.

    """
    obj = model_to_dict(model, sort=sort)
    obj["version"] = YAML_SPEC
    if isinstance(filename, str):
        with io.open(filename, "w") as file_handle:
            yaml.dump(obj, file_handle, **kwargs)
    else:
        yaml.dump(obj, filename, **kwargs)


def load_yaml_model(filename: Union[str, Path]) -> "Model":
    """Load a cobra model from a file in YAML format.

    Parameters
    ----------
    filename : str or file-like
        File path or descriptor that contains the YAML document describing the
        cobra model.

    Returns
    -------
    cobra.Model
        The cobra model as represented in the YAML document.

    See Also
    --------
    from_yaml : Load from a string.

    """
    if isinstance(filename, (str, Path)):
        with io.open(filename, "r") as file_handle:
            return model_from_dict(yaml.load(file_handle))
    else:
        return model_from_dict(yaml.load(filename))
