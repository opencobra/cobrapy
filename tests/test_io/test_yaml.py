"""Test functionalities of I/O in YAML format."""

import json
from pathlib import Path
from typing import Callable

import pytest
from ruamel.yaml import YAML

from cobra import Model
from cobra import io as cio


def test_load_yaml_model(
    compare_models: Callable, cobra_data_directory: Path, mini_model: Model
) -> None:
    """Test the reading of YAML model."""
    yaml_model = cio.load_yaml_model(cobra_data_directory.joinpath("mini.yml"))
    assert compare_models(mini_model, yaml_model) is None


@pytest.mark.xfail(reason="schema outdated")
def test_save_yaml_model(tmp_path: Path, mini_model: Model) -> None:
    """Test the writing of YAML model."""
    jsonschema = pytest.importorskip("jsonschema")
    output_file = tmp_path.joinpath("mini.yml")
    cio.save_yaml_model(mini_model, str(output_file), sort=True)
    # validate against JSONSchema
    yaml = YAML(typ="unsafe")
    with open(output_file, "r") as infile:
        yaml_to_dict = yaml.load(infile)
    dict_to_json = json.dumps(yaml_to_dict)
    loaded = json.loads(dict_to_json)
    assert jsonschema.validate(loaded, cio.json.json_schema)
