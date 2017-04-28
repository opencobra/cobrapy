# -*- coding: utf-8 -*-

from __future__ import absolute_import

import io
import yaml
from collections import OrderedDict

from six import string_types

from cobra.io.dict import model_to_dict, model_from_dict


def ordered_load(stream, Loader=yaml.Loader, object_pairs_hook=OrderedDict):
    """http://stackoverflow.com/a/21912744"""
    class OrderedLoader(Loader):
        pass

    def construct_mapping(loader, node):
        loader.flatten_mapping(node)
        return object_pairs_hook(loader.construct_pairs(node))

    OrderedLoader.add_constructor(
        yaml.resolver.BaseResolver.DEFAULT_MAPPING_TAG, construct_mapping)
    return yaml.load(stream, OrderedLoader)


def ordered_dump(data, stream=None, Dumper=yaml.Dumper, **kwds):
    """http://stackoverflow.com/a/21912744"""
    class OrderedDumper(Dumper):
        pass

    def _dict_representer(dumper, data):
        return dumper.represent_mapping(
            yaml.resolver.BaseResolver.DEFAULT_MAPPING_TAG, data.items())

    OrderedDumper.add_representer(OrderedDict, _dict_representer)
    return yaml.dump(data, stream, OrderedDumper, **kwds)


def to_yaml(model):
    """
    Return the model as a YAML string.
    
    Parameters
    ----------
    model : cobra.Model
    ordered : bool, optional
        Whether to create YAML representation with deterministic order (
        default true).
    """
    return ordered_dump(model_to_dict(model), Dumper=yaml.SafeDumper)


def from_yaml(document):
    """
    Load cobra model from a YAML string
    
    Parameters
    ----------
    document : str
        A string representation of a metabolic model in YAML format.
    """
    return model_from_dict(ordered_load(document, yaml.SafeLoader))


def save_yaml_model(model, filename):
    document = to_yaml(model)
    if isinstance(filename, string_types):
        with io.open(filename, "w") as file_h:
            file_h.write(document)
    else:
        filename.write(document)


def load_yaml_model(filename):
    if isinstance(filename, string_types):
        with io.open(filename, "r") as file_h:
            document = from_yaml(file_h.read())
    else:
        document = from_yaml(filename.read())
    return document
