# -*- coding: utf-8 -*-

"""Test JSON input and output."""

from __future__ import absolute_import

from builtins import open  # Python 2 unicode compatibility.
from os.path import join

import cobra.io as cio

import helpers


def test_from_json(data_directory, mini_model):
    """Test reading a model from a JSON string."""
    with open(join(data_directory, "mini.json"), encoding="utf-8") as handle:
        json_model = cio.from_json(handle.read())
    helpers.assert_equal_models(mini_model, json_model)


def test_load_json_model(data_directory, mini_model):
    """Test reading a model from a JSON file."""
    json_model = cio.load_json_model(join(data_directory, "mini.json"))
    helpers.assert_equal_models(mini_model, json_model)


def test_to_json(data_directory, mini_model):
    """Test writing a model to a JSON string."""
    output = cio.to_json(mini_model, sort=True, pretty=True)
    with open(join(data_directory, "mini.json"), encoding="utf-8") as handle:
        expected = handle.read()
    assert output == expected


def test_save_json_model(tmpdir, data_directory, mini_model):
    """Test writing a model to a JSON file."""
    output_file = tmpdir.join("mini.json")
    cio.save_json_model(mini_model, str(output_file), sort=True, pretty=True)
    # Validate the written file.
    with open(str(output_file), encoding="utf-8") as handle:
        output = handle.read()
    with open(join(data_directory, "mini.json"), encoding="utf-8") as handle:
        expected = handle.read()
    assert output == expected
