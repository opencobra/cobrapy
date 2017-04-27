# -*- coding: utf-8 -*-

from __future__ import absolute_import

from cobra.io.dict import model_to_dict, model_from_dict

"""http://stackoverflow.com/a/21912744"""
import yaml
from collections import OrderedDict

def ordered_load(stream, Loader=yaml.Loader, object_pairs_hook=OrderedDict):
    class OrderedLoader(Loader):
        pass

    def construct_mapping(loader, node):
        loader.flatten_mapping(node)
        return object_pairs_hook(loader.construct_pairs(node))

    OrderedLoader.add_constructor(
        yaml.resolver.BaseResolver.DEFAULT_MAPPING_TAG,
        construct_mapping)
    return yaml.load(stream, OrderedLoader)


def ordered_dump(data, stream=None, Dumper=yaml.Dumper, **kwds):
    class OrderedDumper(Dumper):
        pass

    def _dict_representer(dumper, data):
        return dumper.represent_mapping(
            yaml.resolver.BaseResolver.DEFAULT_MAPPING_TAG,
            data.items())

    OrderedDumper.add_representer(OrderedDict, _dict_representer)
    return yaml.dump(data, stream, OrderedDumper, **kwds)


def to_yaml(model):
    raise NotImplementedError()


def from_yaml(document):
    raise NotImplementedError()


def save_yaml_model(model, filename):
    raise NotImplementedError()


def load_yaml_model(filename):
    raise NotImplementedError()
