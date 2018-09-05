# -*- coding: utf-8 -*-

"""Test functionalities of json.py"""

from __future__ import absolute_import

import json
from os.path import join

import cobra.io as cio
import pytest
from cobra.test.test_io.conftest import compare_models

jsonschema = pytest.importorskip("jsonschema")


@pytest.mark.xfail(reason="schema outdated")
@pytest.mark.skipif(jsonschema is None, reason="jsonschema not installed")
def test_validate_json(data_directory):
    """Validate file according to JSON-schema."""
    with open(join(data_directory, "mini.json"), "r") as infile:
        loaded = json.load(infile)
    assert jsonschema.validate(loaded, cio.json.json_schema)


def test_load_json_model(data_directory, mini_model):
    """Test the reading of JSON model."""
    json_model = cio.load_json_model(join(data_directory, "mini.json"))
    assert compare_models(mini_model, json_model) is None


@pytest.mark.xfail(reason="schema outdated")
@pytest.mark.skipif(jsonschema is None, reason="jsonschema not installed")
def test_save_json_model(tmpdir, mini_model):
    """Test the writing of JSON model."""
    output_file = tmpdir.join("mini.json")
    cio.save_json_model(mini_model, output_file.strpath, pretty=True)
    # validate against JSONSchema
    with open(output_file, "r") as infile:
        loaded = json.load(infile)
    assert jsonschema.validate(loaded, cio.json.json_schema)
