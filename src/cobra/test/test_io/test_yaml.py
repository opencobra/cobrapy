# -*- coding: utf-8 -*-

"""Test functionalities provided by yaml.py"""

from __future__ import absolute_import

import json
from os.path import join

import pytest
from ruamel.yaml import YAML

import cobra.io as cio
from cobra.test.test_io.conftest import compare_models


def test_load_yaml_model(data_directory, mini_model):
    """Test the reading of YAML model."""
    yaml_model = cio.load_yaml_model(join(data_directory, "mini.yml"))
    assert compare_models(mini_model, yaml_model) is None


@pytest.mark.xfail(reason="schema outdated")
def test_save_yaml_model(tmpdir, mini_model):
    jsonschema = pytest.importorskip("jsonschema")
    """Test the writing of YAML model."""
    output_file = tmpdir.join("mini.yml")
    cio.save_yaml_model(mini_model, output_file.strpath, sort=True)
    # validate against schema
    yaml = YAML(typ="unsafe")
    with open(output_file.strpath, "r") as infile:
        yaml_to_dict = yaml.load(infile)
    dict_to_json = json.dumps(yaml_to_dict)
    loaded = json.loads(dict_to_json)
    assert jsonschema.validate(loaded, cio.json.json_schema)
