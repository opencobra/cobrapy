# -*- coding: utf-8 -*-

from __future__ import absolute_import

try:
    import simplejson as json
except ImportError:
    import json
from six import string_types

from cobra.io.dict import model_to_dict, model_from_dict

JSON_SPEC = "1"


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
    filename : str or file-like
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

    if isinstance(filename, string_types):
        with open(filename, "w") as file_handle:
            json.dump(obj, file_handle, **dump_opts)
    else:
        json.dump(obj, filename, **dump_opts)


def load_json_model(filename):
    """
    Load a cobra model from a file in JSON format.

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
    from_json : Load from a string.
    """
    if isinstance(filename, string_types):
        with open(filename, "r") as file_handle:
            return model_from_dict(json.load(file_handle))
    else:
        return model_from_dict(json.load(filename))


json_schema = {
    "$schema": "http://json-schema.org/draft-04/schema#",
    "title": "COBRA",
    "description": "JSON representation of COBRA model",
    "type": "object",
    "properties": {
        "id": {"type": "string"},
        "name": {"type": "string"},
        "description": {"type": "string"},
        "version": {
            "type": "integer",
            "default": 1,
        },

        "reactions": {
            "type": "array",
            "items": {
                "type": "object",
                "properties": {
                    "id": {"type": "string"},
                    "name": {"type": "string"},
                    "metabolites": {
                        "type": "object",
                        "patternProperties": {
                            ".*": {"type": "number"},
                        }
                    },
                    "gene_reaction_rule": {"type": "string"},
                    "lower_bound": {"type": "number"},
                    "upper_bound": {"type": "number"},
                    "objective_coefficient": {
                        "type": "number",
                        "default": 0,
                    },
                    "variable_kind": {
                        "type": "string",
                        "pattern": "integer|continuous",
                        "default": "continuous"
                    },
                    "subsystem": {"type": "string"},
                    "notes": {"type": "object"},
                    "annotation": {"type": "object"},
                },
                "required": ["id", "name", "metabolites", "lower_bound",
                             "upper_bound", "gene_reaction_rule"],
                "additionalProperties": False,
            }
        },
        "metabolites": {
            "type": "array",
            "items": {
                "type": "object",
                "properties": {
                    "id": {"type": "string"},
                    "name": {"type": "string"},
                    "compartment": {
                        "type": "string",
                        "pattern": "[a-z]{1,2}"
                    },
                    "charge": {"type": "integer"},
                    "formula": {"type": "string"},
                    "_bound": {
                        "type": "number",
                        "default": 0
                    },
                    "_constraint_sense": {
                        "type": "string",
                        "default": "E",
                        "pattern": "E|L|G",
                    },
                    "notes": {"type": "object"},
                    "annotation": {"type": "object"},
                },
                "required": ["id", "name", "compartment"],
                "additionalProperties": False,
            }

        },
        "genes": {
            "type": "array",
            "items": {
                "type": "object",
                "properties": {
                    "id": {"type": "string"},
                    "name": {"type": "string"},
                    "notes": {"type": "object"},
                    "annotation": {"type": "object"},
                },
                "required": ["id", "name"],
                "additionalProperties": False,
            }

        },
        "compartments": {
            "type": "object",
            "patternProperties": {
                "[a-z]{1,2}": {"type": "string"}
            }
        },
        "notes": {"type": "object"},
        "annotation": {"type": "object"},
    },
    "required": ["id", "reactions", "metabolites", "genes"],
    "additionalProperties": False,
}
