# -*- coding: utf-8 -*-

"""Test functionalities of json.py"""

from __future__ import absolute_import

from filecmp import cmp
from os.path import join
from warnings import warn

from cobra import io
from cobra.test.test_io.conftest import compare_models

try:
    import jsonschema
except ImportError:
    jsonschema = None


def validate_json(data_directory):
    """Validate file according to JSON-schema."""
    with open(join(data_directory, "mini.json"), "r") as infile:
        loaded = io.json.json.load(infile)
    if jsonschema is None:
        warn("jsonschema not installed")
    else:
        jsonschema.validate(loaded, io.json.json_schema)


def test_load_json_model(data_directory, mini_model):
    """Test the reading of JSON model."""
    json_model = io.load_json_model(join(data_directory, "mini.json"))
    assert compare_models(mini_model, json_model) is None


def test_save_json_model(data_directory, tmpdir):
    """Test the writing of JSON model."""
    input_file = join(data_directory, "mini.json")
    output_file = join(tmpdir, "mini.json")

    json_model = io.load_json_model(input_file)
    io.save_json_model(json_model, output_file, pretty=True)
    assert cmp(input_file, output_file, shallow=False)
