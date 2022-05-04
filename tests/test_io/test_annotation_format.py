"""Test model annotations in JSON format."""

from pathlib import Path

import pytest

from cobra.io import load_json_model, write_sbml_model


def test_load_json_model_valid(data_directory: Path, tmp_path: Path) -> None:
    """Test loading a valid annotation from JSON.

    data_directory : pathlib.Path
        The path to the test data directory.
    tmp_path : pathlib.Path
        The path to the temporary test assets store.

    """
    path_to_file = data_directory / "valid_annotation_format.json"
    model = load_json_model(path_to_file)
    expected = {
        "bigg.reaction": [["is", "PFK26"]],
        "kegg.reaction": [["is", "R02732"]],
        "rhea": [["is", "15656"]],
    }
    for metabolite in model.metabolites:
        assert metabolite.annotation == expected
    path_to_output = tmp_path / "valid_annotation_output.xml"
    write_sbml_model(model, str(path_to_output.resolve()))


def test_load_json_model_invalid(data_directory: Path) -> None:
    """Test that loading an invalid annotation from JSON raises TypeError.

    data_directory : pathlib.Path
        The path to the test data directory.

    """
    path = data_directory / "invalid_annotation_format.json"
    with pytest.raises(TypeError):
        load_json_model(path)
