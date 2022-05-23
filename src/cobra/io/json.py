"""Provide functions for I/O in JSON format."""
from pathlib import Path
from typing import TYPE_CHECKING, Any, Dict, List, Union, IO

from .dict import model_from_dict, model_to_dict
from cobra import io as cio

import jsonschema
from importlib_resources import open_text


try:
    import simplejson as json
except ImportError:
    import json

if TYPE_CHECKING:
    from cobra import Model


JSON_SPEC = "1"


def json_schema_v1() -> Dict:
    with open_text(cio, "schema_v1.json") as handle:
        schema_v1 = json.load(handle)
    return schema_v1


def json_schema_v2() -> Dict:
    with open_text(cio, "schema_v2.json") as handle:
        schema_v2 = json.load(handle)
    return schema_v2


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
    filename: Union[str, Path, IO[str]],
    sort: bool = False,
    pretty: bool = False,
    **kwargs: Any
) -> None:
    """Write the cobra model to a file in JSON format.

    Parameters
    ----------
    model : cobra.Model
        The cobra model to represent.
    filename : str or file-like (IO[str]) or Path
        File path or file handle or str that the JSON representation should be
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


def load_json_model(filename: Union[str, Path, IO[str]]) -> "Model":
    """Load a cobra model from a file in JSON format.

    Parameters
    ----------
    filename : str or file-like or Path
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


def validate_json_model(
    filename: Union[str, Path, IO[str]], json_schema_version: int = 1
) -> List:
    """
    Validate a model in json format against the schema with given version
    Parameters
    ----------
    filename : str or Path or file-like
        File path or descriptor that contains the JSON document describing the
        cobra model.
    json_schema_version : int {1, 2}
        the version of schema to be used for validation.
        Currently we have v1 and v2 only and v2 is under development
    Returns
    -------
    errors : list
        The list of errors encountered while validating
    """

    if json_schema_version == 1:
        schema = json_schema_v1()
    elif json_schema_version == 2:
        schema = json_schema_v2()
    else:
        raise ValueError(
            f"Only v1 and v2 of JSON schema are available. JSON "
            f"schema v{json_schema_version} is not supported."
        )

    # TODO - Should the validator be picked by schema?
    #  Something like validators.validator_for
    validator = jsonschema.Draft7Validator(schema)

    try:
        if isinstance(filename, (str , Path)):
            with open(filename, "r") as file_handle:
                errors = validator.iter_errors(json.load(file_handle))
        else:
            errors = validator.iter_errors(json.load(filename))
    except OSError:
        errors = validator.iter_errors(json.loads(filename))

    return list(errors)
