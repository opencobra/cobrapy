from os.path import join

from cobra.io import load_json_model, write_sbml_model


def test_load_json_model_valid(data_directory, tmp_path):
    """Test loading a valid annotation from JSON."""
    path_to_file = join(data_directory, "valid_annotation_format.json")
    model = load_json_model(path_to_file)
    expected = {
        "bigg.reaction": ["PFK26"],
        "kegg.reaction": ["R02732"],
        "rhea": ["15656"],
    }
    for metabolite in model.metabolites:
        assert metabolite.annotation == expected
    path_to_output = join(str(tmp_path), "valid_annotation_output.xml")
    write_sbml_model(model, path_to_output)


def test_load_json_model_invalid(data_directory):
    """Test that loading an invalid annotation from JSON raises TypeError"""
    path = join(data_directory, "invalid_annotation_format.json")
    # with pytest.raises(TypeError):
    #     model = load_json_model(path)

    # the issue of reading annotation when it is in the form of
    # of list of list has been resolved. When such type of annotation
    # are encountered, they will be first fixed and then added
    model = load_json_model(path)
    anno = model.metabolites[0].annotation
    assert anno == {"kegg.compound": ["C01468"], "chebi": ["CHEBI:11981"]}
