from os.path import join

import pytest

from cobra.core.metadata import Notes
from cobra.io import (
    load_json_model, read_sbml_model, save_json_model, write_sbml_model)


incoming_notes_str = '\
<notes>\n\
  <body xmlns="http://www.w3.org/1999/xhtml">\n\
    <div style="height: 60px; background-color: #09F;">\n\
      <p> Key1 : Value1 </p>\n\
      <p> Key2 : Value2 </p>\n\
      <div style="margin-left: auto; margin-right: auto; ' \
                     'width: 970px;">\n\
        <h1> A Heading </h1>\n\
        <div class="dc:title"> e_coli_core - Escherichia coli ' \
                     'str. K-12 substr. MG1655 </div>\n\
      </div>\n\
      <p> Key3 : Value3 </p>\n\
    </div>\n\
  </body>\n\
</notes>'

modified_notes_str = '\
<notes>\n\
  <body xmlns="http://www.w3.org/1999/xhtml">\n\
    <div style="height: 60px; background-color: #09F;">\n\
      <p> Key1 : New Value 1 </p>\n\
      <p> Key2 : Value2 </p>\n\
      <div style="margin-left: auto; margin-right: auto; ' \
                     'width: 970px;">\n\
        <h1> A Heading </h1>\n\
        <div class="dc:title"> e_coli_core - Escherichia coli ' \
                     'str. K-12 substr. MG1655 </div>\n\
      </div>\n\
      <p> Key3 : New Value 3 </p>\n\
    </div>\n\
  </body>\n\
</notes>'


def test_notes(data_directory, tmp_path):
    """reading notes from SBML to cobra model"""
    model = read_sbml_model(join(data_directory,
                                 "e_coli_core_for_annotation.xml"))
    rx1 = model.reactions[0]
    # making notes object to test equality check of
    # two notes object
    notes = Notes(incoming_notes_str)

    assert rx1.notes.notes_text == incoming_notes_str
    assert rx1.notes == notes

    # keys inside notes dict
    list_of_keys = ['Key1', 'Key2', 'Key3']

    for key in list_of_keys:
        assert key in rx1.notes

    assert rx1.notes['Key1'] == 'Value1'
    assert rx1.notes['Key2'] == 'Value2'
    assert rx1.notes['Key3'] == 'Value3'

    # modifying already present key-value
    rx1.notes["Key1"] = "New Value 1"
    rx1.notes["Key3"] = "New Value 3"

    # trying to insert a new key-value
    with pytest.raises(ValueError):
        rx1.notes["Key4"] = "New Value 3"

    # checking modified notes dict and string
    assert rx1.notes.notes_text == modified_notes_str
    assert rx1.notes['Key1'] == 'New Value 1'
    assert rx1.notes['Key2'] == 'Value2'
    assert rx1.notes['Key3'] == 'New Value 3'

    # writing and reading back the model
    path_to_file = join(str(tmp_path), "model_notes.xml")
    write_sbml_model(model, path_to_file)

    model_after_reading = read_sbml_model(path_to_file)
    rx1_after_reading = model_after_reading.reactions[0]

    # checks after reading model back again
    assert rx1_after_reading.notes.notes_text == modified_notes_str
    assert rx1_after_reading.notes['Key1'] == 'New Value 1'
    assert rx1_after_reading.notes['Key2'] == 'Value2'
    assert rx1_after_reading.notes['Key3'] == 'New Value 3'


def test_reading_writing_notes(data_directory, tmp_path):
    # reading model with notes
    model = read_sbml_model(join(data_directory,
                                 "e_coli_core_for_annotation.xml"))

    # checking notes data
    rx1 = model.reactions[0]
    assert rx1.notes.notes_text == incoming_notes_str

    # reading and writing in json format
    path_to_json = join(str(tmp_path), "json_notes.json")
    save_json_model(model, path_to_json)
    model_from_json = load_json_model(path_to_json)
    rx1_from_json = model_from_json.reactions[0]
    assert rx1_from_json.notes.notes_text == incoming_notes_str


# def test_notes(tmp_path):
#     """Testing if model notes are written in SBML"""
#     path_to_file = join(str(tmp_path), "model_notes.xml")

#     # making a minimal cobra model to test notes
#     model = cobra.Model("e_coli_core")
#     model.notes["Remark"] = "...Model Notes..."
#     met = cobra.Metabolite("pyr_c", compartment="c")
#     model.add_metabolites([met])
#     met.notes["Remark"] = "Note with \n newline"
#     rxn = cobra.Reaction("R_ATPM")
#     model.add_reactions([rxn])
#     rxn.notes["Remark"] = "What about me?"
#     model.objective_direction = "max"
#     model.objective = rxn
#     write_sbml_model(model, path_to_file)

#     # reading the model back
#     model_after_reading = read_sbml_model(path_to_file)
#     met_after_reading = model_after_reading.metabolites.get_by_id("pyr_c")
#     reaction_after_reading =
#     model_after_reading.reactions.get_by_id("R_ATPM")

#     # checking if notes are written to model
#     assert model_after_reading.notes["Remark"] == "...Model Notes..."

#     # checking notes for metabolite and reaction
#     assert met_after_reading.notes["Remark"] == "Note with \n newline"
#     assert reaction_after_reading.notes["Remark"] == "What about me?"
