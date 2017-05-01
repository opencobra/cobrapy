# -*- coding: utf-8 -*-

from __future__ import absolute_import

try:
    import simplejson as json
except ImportError:
    import json
from six import string_types

from cobra.io.dict import model_to_dict, model_from_dict

JSON_SPEC = "1"


def to_json(model):
    """Save the cobra model as a json string"""
    obj = model_to_dict(model)
    obj["version"] = JSON_SPEC
    return json.dumps(obj, allow_nan=False)


def from_json(jsons):
    """Load cobra model from a json string"""
    return model_from_dict(json.loads(jsons))


def load_json_model(file_name):
    """Load a cobra model stored as a json file

    Parameters
    ----------
    file_name : str or file-like object

    Returns
    -------
    cobra.Model
       The loaded model
    """
    # open the file
    should_close = False
    if isinstance(file_name, string_types):
        file_name = open(file_name, 'r')
        should_close = True

    model = model_from_dict(json.load(file_name))

    if should_close:
        file_name.close()

    return model


def save_json_model(model, file_name, pretty=False):
    """Save the cobra model as a json file.

    Parameters
    ----------
    model : cobra.core.Model.Model
        The model to save
    file_name : str or file-like object
        The file to save to
    """
    # open the file
    obj = model_to_dict(model)
    obj["version"] = JSON_SPEC
    should_close = False
    if isinstance(file_name, string_types):
        file_name = open(file_name, 'w')
        should_close = True

    if pretty:
        dump_opts = {"indent": 4, "separators": (",", ": "), "sort_keys": True}
    else:
        dump_opts = {}

    json.dump(obj, file_name, allow_nan=False, **dump_opts)

    if should_close:
        file_name.close()


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
