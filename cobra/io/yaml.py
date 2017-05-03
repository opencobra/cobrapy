# -*- coding: utf-8 -*-

from __future__ import absolute_import

import io

from six import string_types
from ruamel import yaml

from cobra.io.dict import model_to_dict, model_from_dict

YAML_SPEC = "1"


def to_yaml(model, **kwargs):
    """
    Return the model as a YAML document.

    ``kwargs`` are passed on to ``yaml.dump``.

    Parameters
    ----------
    model : cobra.Model
        The cobra model to represent.

    Returns
    -------
    str
        String representation of the cobra model as a YAML document.

    See Also
    --------
    save_yaml_model : Write directly to a file.
    ruamel.yaml.dump : Base function.
    """
    obj = model_to_dict(model)
    obj["version"] = YAML_SPEC
    return yaml.dump(obj, Dumper=yaml.RoundTripDumper, **kwargs)


def from_yaml(document):
    """
    Load a cobra model from a YAML document.

    Parameters
    ----------
    document : str
        The YAML document representation of a cobra model.

    Returns
    -------
    cobra.Model
        The cobra model as represented in the YAML document.

    See Also
    --------
    load_yaml_model : Load directly from a file.
    """
    return model_from_dict(yaml.load(document, yaml.RoundTripLoader))


def save_yaml_model(model, filename, **kwargs):
    """
    Write the cobra model to a file in YAML format.

    ``kwargs`` are passed on to ``yaml.dump``.

    Parameters
    ----------
    model : cobra.Model
        The cobra model to represent.
    filename : str or file-like
        File path or descriptor that the YAML representation should be
        written to.

    See Also
    --------
    to_yaml : Return a string representation.
    ruamel.yaml.dump : Base function.
    """
    obj = model_to_dict(model)
    obj["version"] = YAML_SPEC
    if isinstance(filename, string_types):
        with io.open(filename, "w") as file_handle:
            yaml.dump(obj, file_handle, Dumper=yaml.RoundTripDumper, **kwargs)
    else:
        yaml.dump(obj, filename, Dumper=yaml.RoundTripDumper, **kwargs)


def load_yaml_model(filename):
    """
    Load a cobra model from a file in YAML format.

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
    if isinstance(filename, string_types):
        with io.open(filename, "r") as file_handle:
            return model_from_dict(yaml.load(file_handle,
                                             yaml.RoundTripLoader))
    else:
        return model_from_dict(yaml.load(filename, yaml.RoundTripLoader))
