# -*- coding: utf-8 -*-

"""Test functionalities of json.py"""

from os.path import join

from cobra import io as cio
from cobra.test.test_io.conftest import compare_models


def test_validate_json(data_directory):
    """Validate file according to JSON-schema."""
    path_old_format = join(data_directory, "e_coli_core.json")
    # validate the model using JSON schema v1
    assert cio.validate_json_model(filename=path_old_format,
                                   json_schema_version=1) == (True, "")

    path_new_format = join(data_directory, "e_coli_new_format.json")
    # validate the model using JSON schema v2
    assert cio.validate_json_model(filename=path_new_format,
                                   json_schema_version=2) == (True, "")


def test_load_json_model(data_directory, mini_model):
    """Test the reading of JSON model."""
    json_model = cio.load_json_model(join(data_directory, "mini.json"))
    assert compare_models(mini_model, json_model) is None


def test_save_json_model(tmpdir, mini_model):
    """Test the writing of JSON model."""
    output_file = tmpdir.join("mini.json")
    cio.save_json_model(mini_model, output_file.strpath, pretty=True)
    # validate against JSONSchema
    assert cio.validate_json_model(output_file, 1) == (True, "")


def test_reaction_bounds_json(data_directory, tmp_path):
    """Test reading and writing of model with inf bounds in json"""

    """Path to XML file with INF bounds"""
    path_to_xml_inf_file = join(data_directory, "fbc_ex1.xml")
    model_xml_inf = cio.read_sbml_model(path_to_xml_inf_file)
    path_to_output = join(str(tmp_path), "fbc_ex1_json.json")

    """Saving model with inf bounds in json form without error"""
    cio.save_json_model(model_xml_inf, path_to_output)

    """Path to JSON file with INF bounds"""
    path_to_JSON_inf_file = join(data_directory, "JSON_with_inf_bounds.json")
    model_json_inf = cio.load_json_model(path_to_JSON_inf_file)
    assert model_json_inf.reactions[0].upper_bound == float("inf")
