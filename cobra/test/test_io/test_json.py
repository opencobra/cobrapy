# -*- coding: utf-8 -*-

"""Test functionalities of json.py"""

from __future__ import absolute_import

import json
from os.path import join

import pytest

import cobra.io as cio
from cobra.test.test_io.conftest import compare_models


def test_validate_json(data_directory, json_schema_v1):
    """Validate file according to JSON-schema."""
    jsonschema = pytest.importorskip("jsonschema")
    with open(join(data_directory, "mini.json"),
              "r") as infile:
        loaded = json.load(infile)
    assert jsonschema.validate(loaded, json_schema_v1) is None


def test_load_json_model(data_directory, mini_model):
    """Test the reading of JSON model."""
    json_model = cio.load_json_model(join(data_directory, "mini.json"))
    assert compare_models(mini_model, json_model) is None


def test_save_json_model(tmpdir, mini_model, json_schema_v1):
    """Test the writing of JSON model."""
    jsonschema = pytest.importorskip("jsonschema")
    output_file = tmpdir.join("mini.json")
    cio.save_json_model(mini_model, output_file.strpath, pretty=True)
    # validate against JSONSchema
    with open(str(output_file), "r") as infile:
        loaded = json.load(infile)
    assert jsonschema.validate(loaded, json_schema_v1) is None
