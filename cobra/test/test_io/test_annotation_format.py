from os.path import join

import pytest

from cobra.io import load_json_model, write_sbml_model


def test_load_json_model_valid(data_directory, tmp_path):
    """Test loading a valid annotation from JSON."""
    path_to_file = join(data_directory, "valid_annotation_format.json")
    model = load_json_model(path_to_file)
    expected = {
        'bigg.reaction': [['is', 'PFK26']],
        'kegg.reaction': [['is', 'R02732']],
        'rhea': [['is', '15656']]
    }
    for metabolite in model.metabolites:
        assert metabolite.annotation == expected
    path_to_output = join(str(tmp_path), 'valid_annotation_output.xml')
    write_sbml_model(model, path_to_output)


def test_load_json_model_invalid(data_directory):
    """Test that loading an invalid annotation from JSON raises TypeError"""
    path = join(data_directory, "invalid_annotation_format.json")
    with pytest.raises(TypeError):
        model = load_json_model(path)
