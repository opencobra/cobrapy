# -*- coding: utf-8 -*-

"""Test functionalities provided by yaml.py"""

import json
from os.path import join

from ruamel.yaml import YAML

from cobra import io as cio
from cobra.test.test_io.conftest import compare_models


def test_load_yaml_model(data_directory, mini_model):
    """Test the reading of YAML model."""
    yaml_model = cio.load_yaml_model(join(data_directory, "mini.yml"))
    assert compare_models(mini_model, yaml_model) is None


def test_save_yaml_model(tmpdir, mini_model):
    """Test the writing of YAML model."""
    output_file = tmpdir.join("mini.yml")
    cio.save_yaml_model(mini_model, output_file.strpath, sort=True)
    # validate against schema
    yaml = YAML(typ="unsafe")
    with open(output_file.strpath, "r") as infile:
        yaml_to_dict = yaml.load(infile)
    dict_to_json = json.dumps(yaml_to_dict)
    errors = cio.validate_json_model(filename=dict_to_json, json_schema_version=1)
    assert len(errors) == 0
