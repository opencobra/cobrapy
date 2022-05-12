"""Test functionalities of I/O in JSON format."""

from os.path import join
from pathlib import Path
from typing import Callable


from cobra import Model
from cobra import io as cio


def test_validate_json(data_directory: Path) -> None:
    """Validate file according to JSON-schema."""
    path_old_format = join(data_directory, "e_coli_core.json")
    # validate the model using JSON schema v1
    list_errors = cio.validate_json_model(
        filename=path_old_format, json_schema_version=1
    )
    assert len(list_errors) == 0

    path_new_format = join(data_directory, "e_coli_new_format.json")
    # validate the model using JSON schema v2
    errors = cio.validate_json_model(filename=path_new_format, json_schema_version=2)
    assert len(errors) == 0

    # test for invalid json model according to schema
    errors_invalid = cio.validate_json_model(
        filename=path_old_format, json_schema_version=2
    )
    assert len(errors_invalid) == 309


def test_load_json_model(
    compare_models: Callable, data_directory: Path, mini_model: Model
) -> None:
    """Test the reading of JSON model."""
    json_model = cio.load_json_model(data_directory / "mini.json")
    assert compare_models(mini_model, json_model) is None


def test_save_json_model(
    tmp_path: Path,
    mini_model: Model
) -> None:
    """Test the writing of JSON model."""
    output_file = tmp_path.joinpath("mini.json")
    cio.save_json_model(mini_model, output_file, pretty=True)
    # validate against JSONSchema
    errors = cio.validate_json_model(output_file, 1)
    assert len(errors) == 0


def test_reaction_bounds_json(data_directory: Path, tmp_path: Path) -> None:
    """Test reading and writing of model with inf bounds in JSON."""
    """Path to XML file with INF bounds"""
    path_to_xml_inf_file = join(data_directory, "fbc_ex1.xml")
    model_xml_inf = cio.read_sbml_model(path_to_xml_inf_file)
    path_to_output = tmp_path.joinpath("fbc_ex1_json.json")

    """Saving model with inf bounds in json form without error"""
    cio.save_json_model(model_xml_inf, path_to_output)

    """Path to JSON file with INF bounds"""
    path_to_JSON_inf_file = data_directory.joinpath("JSON_with_inf_bounds.json")
    model_json_inf = cio.load_json_model(path_to_JSON_inf_file)
    assert model_json_inf.reactions[0].upper_bound == float("inf")
