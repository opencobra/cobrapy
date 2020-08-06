# -*- coding: utf-8 -*-

from pathlib import Path

import jsonschema
from importlib_resources import open_text
from six import string_types

from cobra import io as cio
from cobra.io.dict import model_from_dict, model_to_dict


try:
    import simplejson as json
except ImportError:
    import json


JSON_SPEC = "1"


def json_schema_v1():
    with open_text(cio, "schema_v1.json") as handle:
        schema_v1 = json.load(handle)
    return schema_v1


def json_schema_v2():
    with open_text(cio, "schema_v2.json") as handle:
        schema_v2 = json.load(handle)
    return schema_v2


def to_json(model, sort=False, **kwargs):
    """
    Return the model as a JSON document.

    ``kwargs`` are passed on to ``json.dumps``.

    Parameters
    ----------
    model : cobra.Model
        The cobra model to represent.
    sort : bool, optional
        Whether to sort the metabolites, reactions, and genes or maintain the
        order defined in the model.

    Returns
    -------
    str
        String representation of the cobra model as a JSON document.

    See Also
    --------
    save_json_model : Write directly to a file.
    json.dumps : Base function.
    """
    obj = model_to_dict(model, sort=sort)
    obj[u"version"] = JSON_SPEC
    return json.dumps(obj, allow_nan=False, **kwargs)


def from_json(document):
    """
    Load a cobra model from a JSON document.

    Parameters
    ----------
    document : str
        The JSON document representation of a cobra model.

    Returns
    -------
    cobra.Model
        The cobra model as represented in the JSON document.

    See Also
    --------
    load_json_model : Load directly from a file.
    """
    return model_from_dict(json.loads(document))


def save_json_model(model, filename, sort=False, pretty=False, **kwargs):
    """
    Write the cobra model to a file in JSON format.

    ``kwargs`` are passed on to ``json.dump``.

    Parameters
    ----------
    model : cobra.Model
        The cobra model to represent.
    filename : str or file-like or Path
        File path or descriptor that the JSON representation should be
        written to.
    sort : bool, optional
        Whether to sort the metabolites, reactions, and genes or maintain the
        order defined in the model.
    pretty : bool, optional
        Whether to format the JSON more compactly (default) or in a more
        verbose but easier to read fashion. Can be partially overwritten by the
        ``kwargs``.

    See Also
    --------
    to_json : Return a string representation.
    json.dump : Base function.
    """
    obj = model_to_dict(model, sort=sort)
    obj[u"version"] = JSON_SPEC

    if pretty:
        dump_opts = {
            "indent": 4, "separators": (",", ": "), "sort_keys": True,
            "allow_nan": False}
    else:
        dump_opts = {
            "indent": 0, "separators": (",", ":"), "sort_keys": False,
            "allow_nan": False}
    dump_opts.update(**kwargs)

    if isinstance(filename, (string_types, Path)):
        with open(filename, "w") as file_handle:
            json.dump(obj, file_handle, **dump_opts)
    else:
        json.dump(obj, filename, **dump_opts)


def load_json_model(filename):
    """
    Load a cobra model from a file in JSON format.

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
    from_json : Load from a string.
    """
    if isinstance(filename, (string_types, Path)):
        with open(filename, "r") as file_handle:
            return model_from_dict(json.load(file_handle))
    else:
        return model_from_dict(json.load(filename))


def validate_json_model(filename, json_schema_version=1):
    """
    Validate a model in json format against the schema with given version
    Parameters
    ----------
    filename : str or file-like
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

    if json_schema_version not in [1, 2]:
        return False, "Incorrect version passed for JSON schema. COBRApy \
                       only supports v1 and v2 of JSON schema"

    if json_schema_version == 1:
        schema = json_schema_v1()
    elif json_schema_version == 2:
        schema = json_schema_v2()
    else:
        raise ValueError("Only v1 and v2 of JSON schema are available. JSON "
                         "schema v{} is not supported.".format(json_schema_version))

    validator = jsonschema.Draft7Validator(schema)

    if isinstance(filename, string_types):
        with open(filename, "r") as file_handle:
            errors = validator.iter_errors(json.load(file_handle))
    else:
        errors = validator.iter_errors(json.load(filename))

    return list(errors)
