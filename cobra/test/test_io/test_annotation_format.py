import pytest
from os.path import join
from cobra.io import load_json_model, write_sbml_model


def test_load_json_model_valid(data_directory):
    """Testing valid annotation format"""
    path_to_file = join(data_directory, "valid_annotation_format.json")
    model = load_json_model(path_to_file)
    dict = {
        'bigg.reaction': [['is', 'PFK26']],
        'kegg.reaction': [['is', 'R02732']],
        'rhea': [['is', '15656']]
    }
    for x in model.metabolites:
        assert x.annotation == dict, "Dictionary is not equal"
    path_to_output = join(data_directory, 'valid_annotation_output.xml')
    write_sbml_model(model, path_to_output)


def test_load_json_model_invalid(data_directory):
    """Testing invalid annotation format"""
    path = join(data_directory, "invalid_annotation_format.json")
    with pytest.raises(TypeError):
        model = load_json_model(path)
