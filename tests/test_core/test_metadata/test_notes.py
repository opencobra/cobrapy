import os
from pathlib import Path

import pytest

from cobra import Metabolite, Model, Reaction
from cobra.core import Notes
from cobra.io import load_json_model, read_sbml_model, save_json_model, write_sbml_model


def test_notes_io(tmp_path: Path) -> None:
    """Test if model notes are written and read from/to SBML.

    Parameters
    ----------
    tmp_path : pathlib.Path
        The path to the temporary test assets store.

    """
    path_to_file = tmp_path / "model_notes.xml"

    # making a minimal cobra model to test notes
    model = Model("e_coli_core")
    model.notes = Notes().notes_from_dict({"Remark": "...Model Notes..."})
    met = Metabolite("pyr_c", compartment="c")
    model.add_metabolites([met])
    met.notes = Notes().notes_from_dict({"Remark": "Note with \n newline"})
    rxn = Reaction("R_ATPM")
    model.add_reactions([rxn])
    rxn.notes = Notes().notes_from_dict({"Remark": "What about me?"})
    model.objective_direction = "max"
    model.objective = rxn
    write_sbml_model(model, str(path_to_file.resolve()))

    # reading the model back
    model_after_reading = read_sbml_model(str(path_to_file.resolve()))
    met_after_reading = model_after_reading.metabolites.get_by_id("pyr_c")
    reaction_after_reading = model_after_reading.reactions.get_by_id("R_ATPM")

    # checking if notes are written to model
    assert model_after_reading.notes["Remark"] == "...Model Notes..."

    # checking notes for metabolite and reaction
    assert met_after_reading.notes["Remark"] == "Note with \n newline"
    assert reaction_after_reading.notes["Remark"] == "What about me?"


NEW_VALUE1 = "New Value 1"
NEW_VALUE3 = "New Value 3"


incoming_notes_str = (
    '\
<notes>\n\
  <body xmlns="http://www.w3.org/1999/xhtml">\n\
    <div style="height: 60px; background-color: #09F;">\n\
      <p> Key1 : Value1 </p>\n\
      <p> Key2 : Value2 </p>\n\
      <div style="margin-left: auto; margin-right: auto; '
    'width: 970px;">\n\
        <h1> A Heading </h1>\n\
        <div class="dc:title"> e_coli_core - Escherichia coli '
    "str. K-12 substr. MG1655 </div>\n\
      </div>\n\
      <p> Key3 : Value3 </p>\n\
    </div>\n\
  </body>\n\
</notes>"
)

modified_notes_str = (
    '\
<notes>\n\
  <body xmlns="http://www.w3.org/1999/xhtml">\n\
    <div style="height: 60px; background-color: #09F;">\n\
      <p> Key1 : New Value 1 </p>\n\
      <p> Key2 : Value2 </p>\n\
      <div style="margin-left: auto; margin-right: auto; '
    'width: 970px;">\n\
        <h1> A Heading </h1>\n\
        <div class="dc:title"> e_coli_core - Escherichia coli '
    "str. K-12 substr. MG1655 </div>\n\
      </div>\n\
      <p> Key3 : New Value 3 </p>\n\
    </div>\n\
  </body>\n\
</notes>"
)


def test_notes(data_directory, tmp_path):
    """reading notes from SBML to cobra model"""
    model_path = os.path.join(data_directory, "e_coli_core_for_annotation.xml")
    assert os.path.exists(model_path)
    model = read_sbml_model(model_path)
    rx1 = model.reactions[0]
    # making notes object to test equality check of
    # two notes object
    notes = Notes(incoming_notes_str)

    assert rx1.notes.notes_xhtml == incoming_notes_str
    assert rx1.notes == notes

    # keys inside notes dict
    list_of_keys = ["Key1", "Key2", "Key3"]

    for key in list_of_keys:
        assert key in rx1.notes

    assert rx1.notes["Key1"] == "Value1"
    assert rx1.notes["Key2"] == "Value2"
    assert rx1.notes["Key3"] == "Value3"

    # modifying already present key-value
    rx1.notes["Key1"] = NEW_VALUE1
    rx1.notes["Key3"] = NEW_VALUE3

    # trying to insert a new key-value
    with pytest.raises(ValueError):
        rx1.notes["Key4"] = NEW_VALUE3

    # checking modified notes dict and string
    assert rx1.notes.notes_xhtml == modified_notes_str
    assert rx1.notes["Key1"] == NEW_VALUE1
    assert rx1.notes["Key2"] == "Value2"
    assert rx1.notes["Key3"] == NEW_VALUE3

    # writing and reading back the model
    path_to_file = os.path.join(tmp_path, "model_notes.xml")
    write_sbml_model(model, path_to_file)

    model_after_reading = read_sbml_model(path_to_file)
    rx1_after_reading = model_after_reading.reactions[0]

    # checks after reading model back again
    assert rx1_after_reading.notes.notes_xhtml == modified_notes_str
    assert rx1_after_reading.notes["Key1"] == NEW_VALUE1
    assert rx1_after_reading.notes["Key2"] == "Value2"
    assert rx1_after_reading.notes["Key3"] == NEW_VALUE3


def test_reading_writing_notes(data_directory, tmp_path):
    # reading model with notes
    model = read_sbml_model(
        os.path.join(data_directory, "e_coli_core_for_annotation.xml")
    )

    # checking notes data
    rx1 = model.reactions[0]
    assert rx1.notes.notes_xhtml == incoming_notes_str

    # reading and writing in json format
    path_to_json = os.path.join(str(tmp_path), "json_notes.json")
    save_json_model(model, path_to_json)
    model_from_json = load_json_model(path_to_json)
    rx1_from_json = model_from_json.reactions[0]
    assert rx1_from_json.notes.notes_xhtml == incoming_notes_str
