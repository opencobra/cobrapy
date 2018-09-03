# -*- coding: utf-8 -*-

"""Test functionalities of json.py"""

from __future__ import absolute_import

from os.path import join
from warnings import warn

import pytest
from cobra import io
from cobra.test.test_io.conftest import compare_models

try:
    import jsonschema
except ImportError:
    jsonschema = None


@pytest.mark.xfail(reason="schema outdated")
def test_validate_json(data_directory):
    """Validate file according to JSON-schema."""
    with open(join(data_directory, "mini.json"), "r") as infile:
        loaded = io.json.json.load(infile)
    if jsonschema is None:
        warn("jsonschema not installed")
    else:
        assert jsonschema.validate(loaded, io.json.json_schema)


def test_load_json_model(data_directory, mini_model):
    """Test the reading of JSON model."""
    json_model = io.load_json_model(join(data_directory, "mini.json"))
    assert compare_models(mini_model, json_model) is None


def test_save_json_model(tmpdir, mini_model):
    """Test the writing of JSON model."""
    output_file = tmpdir.join("mini.json")
    io.save_json_model(mini_model, output_file, pretty=True)
    assert output_file.check()
