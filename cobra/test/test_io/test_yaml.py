# -*- coding: utf-8 -*-

"""Test YAML input and output."""

from __future__ import absolute_import

from builtins import open  # Python 2 unicode compatibility.
from os.path import join

import cobra.io as cio

import helpers


def test_from_yaml(data_directory, mini_model):
    """Test reading a model from a YAML string."""
    with open(join(data_directory, "mini.yml"), encoding="utf-8") as handle:
        yaml_model = cio.from_yaml(handle.read())
    helpers.assert_equal_models(mini_model, yaml_model)


def test_load_yaml_model(data_directory, mini_model):
    """Test reading a model from a YAML file."""
    yaml_model = cio.load_yaml_model(join(data_directory, "mini.yml"))
    helpers.assert_equal_models(mini_model, yaml_model)


def test_to_yaml(data_directory, mini_model):
    """Test writing a model to a YAML string."""
    output = cio.to_yaml(mini_model, sort=True)
    with open(join(data_directory, "mini.yml"), encoding="utf-8") as handle:
        expected = handle.read()
    assert output == expected


def test_save_yaml_model(tmpdir, data_directory, mini_model):
    """Test writing a model to a YAML file."""
    output_file = tmpdir.join("mini.yml")
    cio.save_yaml_model(mini_model, str(output_file), sort=True)
    # Validate the written file.
    with open(str(output_file), encoding="utf-8") as handle:
        output = handle.read()
    with open(join(data_directory, "mini.yml"), encoding="utf-8") as handle:
        expected = handle.read()
    assert output == expected
