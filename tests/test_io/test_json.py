"""Test functionalities of I/O in JSON format."""

import json
from pathlib import Path
from typing import Any, Callable, Dict, Union

import pytest
from importlib_resources import files

from cobra import Model
from cobra import io as cio


@pytest.fixture(scope="module")
def json_schema_v1() -> Dict[str, Union[str, bool, Any]]:
    """Fixture for cobra JSON-schema."""
    with files(cio).joinpath("schema_v1.json").open("r") as handle:
        schema_v1 = json.load(handle)
    return schema_v1


def test_validate_json(
    cobra_data_directory: Path, json_schema_v1: Dict[str, Union[str, bool, Any]]
) -> None:
    """Validate file according to JSON-schema."""
    jsonschema = pytest.importorskip("jsonschema")
    with open(
        cobra_data_directory.joinpath("mini.json"), "r", encoding="utf-8"
    ) as infile:
        loaded = json.load(infile)
    assert jsonschema.validate(loaded, json_schema_v1) is None


def test_load_json_model(
    compare_models: Callable, cobra_data_directory: Path, mini_model: Model
) -> None:
    """Test the reading of JSON model."""
    json_model = cio.load_json_model(cobra_data_directory / "mini.json")
    assert compare_models(mini_model, json_model) is None
    json_model = cio.load_json_model(str(cobra_data_directory / "mini.json"))
    assert compare_models(mini_model, json_model) is None
    with open(cobra_data_directory / "mini.json", "r") as json_handle:
        json_model = cio.load_json_model(json_handle)
        assert compare_models(mini_model, json_model) is None


def test_save_json_model(
    tmp_path: Path,
    mini_model: Model,
    json_schema_v1: Dict[str, Union[str, bool, Any]],
) -> None:
    """Test the writing of JSON model."""
    jsonschema = pytest.importorskip("jsonschema")
    output_file = tmp_path.joinpath("mini.json")
    cio.save_json_model(mini_model, output_file, pretty=True)
    # validate against JSONSchema
    with open(output_file, "r") as infile:
        loaded = json.load(infile)
    assert jsonschema.validate(loaded, json_schema_v1) is None
    output_file.unlink()
    cio.save_json_model(mini_model, str(output_file), pretty=True)
    # validate against JSONSchema
    with open(output_file, "r") as infile:
        loaded = json.load(infile)
    assert jsonschema.validate(loaded, json_schema_v1) is None
    output_file.unlink()
    with output_file.open("w+") as json_outfile:
        cio.save_json_model(mini_model, json_outfile, pretty=True)
        # validate against JSONSchema
        json_outfile.seek(0, 0)
        loaded = json.load(json_outfile)
        assert jsonschema.validate(loaded, json_schema_v1) is None


def test_reaction_bounds_json(data_directory: Path, tmp_path: Path) -> None:
    """Test reading and writing of model with inf bounds in JSON."""
    # Path to XML file with INF bounds
    path_to_xml_inf_file = data_directory / "fbc_ex1.xml"
    model_xml_inf = cio.read_sbml_model(path_to_xml_inf_file)
    path_to_output = tmp_path.joinpath("fbc_ex1_json.json")
    # Saving model with inf bounds in json form without error
    cio.save_json_model(model_xml_inf, path_to_output)
    # Path to JSON file with INF bounds
    path_to_JSON_inf_file = data_directory.joinpath("JSON_with_inf_bounds.json")
    model_json_inf = cio.load_json_model(path_to_JSON_inf_file)
    assert model_json_inf.reactions[0].upper_bound == float("inf")
