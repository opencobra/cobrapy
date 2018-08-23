# -*- coding: utf-8 -*-

"""Test functionalities provided by yaml.py"""

from __future__ import absolute_import

from filecmp import cmp
from os.path import join

from cobra import io
from cobra.test.test_io.conftest import compare_models


def test_load_yaml_model(data_directory, mini_model):
    """Test the reading of YAML model."""
    yaml_model = io.load_yaml_model(join(data_directory, "mini.yml"))
    assert compare_models(mini_model, yaml_model) is None


def test_save_yaml_model(data_directory, tmpdir):
    """Test the writing of YAML model."""
    input_file = join(data_directory, "mini.yml")
    output_file = join(tmpdir, "mini.yml")

    yaml_model = io.load_yaml_model(input_file)
    io.save_yaml_model(yaml_model, output_file, sort=True)
    assert cmp(input_file, output_file, shallow=False)
