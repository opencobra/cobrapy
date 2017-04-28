# -*- coding: utf-8 -*-

from __future__ import absolute_import

import io
from collections import OrderedDict

from six import string_types, iteritems
from ruamel import yaml

from cobra.io.dict import model_to_dict, model_from_dict


def to_yaml(model, **kwargs):
    """
    Return the model as a YAML string.

    Parameters
    ----------
    model : cobra.Model
    """
    return yaml.dump(model_to_dict(model), Dumper=yaml.RoundTripDumper,
                     **kwargs)


def from_yaml(document):
    """
    Load cobra model from a YAML string

    Parameters
    ----------
    document : str
        A string representation of a metabolic model in YAML format.
    """
    return model_from_dict(yaml.load(document, yaml.RoundTripLoader))


def save_yaml_model(model, filename, **kwargs):
    document = to_yaml(model)
    if isinstance(filename, string_types):
        with io.open(filename, "w") as file_h:
            yaml.dump(model_to_dict(model), file_h,
                      Dumper=yaml.RoundTripDumper, **kwargs)
    else:
        yaml.dump(model_to_dict(model), filename, Dumper=yaml.RoundTripDumper,
                  **kwargs)


def load_yaml_model(filename):
    if isinstance(filename, string_types):
        with io.open(filename, "r") as file_h:
            return model_from_dict(yaml.load(file_h, yaml.RoundTripLoader))
    else:
        return model_from_dict(yaml.load(filename, yaml.RoundTripLoader))
