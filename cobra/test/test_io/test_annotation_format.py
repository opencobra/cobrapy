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
    """Test loading an annotation in the form of list of list"""
    path = join(data_directory, "invalid_annotation_format.json")
    expected = {
        'kegg.compound': ['C01468'],
        'chebi': ['CHEBI:11981', 'CHEBI:17847']
    }
    model = load_json_model(path)
    assert model.metabolites[0].annotation == expected
