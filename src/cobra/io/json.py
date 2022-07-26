"""Provide functions for I/O in JSON format."""
from pathlib import Path
from typing import IO, TYPE_CHECKING, Any, Union

from .dict import model_from_dict, model_to_dict


try:
    import simplejson as json
except ImportError:
    import json

if TYPE_CHECKING:
    from cobra import Model


JSON_SPEC = "1"


def to_json(model: "Model", sort: bool = False, **kwargs: Any) -> str:
    """Return the model as a JSON string.

    Parameters
    ----------
    model : cobra.Model
        The cobra model to represent.
    sort : bool, optional
        Whether to sort the metabolites, reactions, and genes or maintain the
        order defined in the model (default False).
    **kwargs : Any
        Keyword arguments passed on to `json.dumps`.

    Returns
    -------
    str
        JSON string representation of the cobra model.

    See Also
    --------
    save_json_model : Write directly to a file.
    json.dumps : Base function.

    """
    obj = model_to_dict(model, sort=sort)
    obj["version"] = JSON_SPEC
    return json.dumps(obj, allow_nan=False, **kwargs)


def from_json(document: str) -> "Model":
    """Load a cobra model from a JSON string.

    Parameters
    ----------
    document : str
        The JSON string representation of a cobra model.

    Returns
    -------
    cobra.Model
        The cobra model as interpreted from the JSON string.

    See Also
    --------
    load_json_model : Load directly from a file.

    """
    return model_from_dict(json.loads(document))


def save_json_model(
    model: "Model",
    filename: Union[str, Path, IO],
    sort: bool = False,
    pretty: bool = False,
    **kwargs: Any
) -> None:
    """Write the cobra model to a file in JSON format.

    Parameters
    ----------
    model : cobra.Model
        The cobra model to represent.
    filename : str or file-like
        File path or descriptor that the JSON representation should be
        written to.
    sort : bool, optional
        Whether to sort the metabolites, reactions, and genes or maintain the
        order defined in the model (default False).
    pretty : bool, optional
        Whether to format the JSON more compactly (default) or in a more
        verbose but easier to read fashion. Can be partially overwritten by the
        `**kwargs` (default False).
    **kwargs : Any
        Keyword arguments passed to `json.dump`.

    See Also
    --------
    to_json : Return a string representation.
    json.dump : Base function.

    """
    obj = model_to_dict(model, sort=sort)
    obj["version"] = JSON_SPEC

    if pretty:
        dump_opts = {
            "indent": 4,
            "separators": (",", ": "),
            "sort_keys": True,
            "allow_nan": False,
        }
    else:
        dump_opts = {
            "indent": 0,
            "separators": (",", ":"),
            "sort_keys": False,
            "allow_nan": False,
        }
    dump_opts.update(**kwargs)

    if isinstance(filename, (str, Path)):
        with open(filename, "w") as file_handle:
            json.dump(obj, file_handle, **dump_opts)
    else:
        json.dump(obj, filename, **dump_opts)


def load_json_model(filename: Union[str, Path, IO]) -> "Model":
    """Load a cobra model from a file in JSON format.

    Parameters
    ----------
    filename : str or file-like
        File path or descriptor that contains the JSON document describing the
        cobra model.

    Returns
    -------
    cobra.Model
        The cobra model as represented in the JSON document.

    See Also
    --------
    from_json : Load from a JSON string.

    """
    if isinstance(filename, (str, Path)):
        with open(filename, "r") as file_handle:
            return model_from_dict(json.load(file_handle))
    else:
        return model_from_dict(json.load(filename))
