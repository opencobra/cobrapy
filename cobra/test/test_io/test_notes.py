from os.path import join

import pytest

import cobra
from cobra.io import read_sbml_model, write_sbml_model


def test_notes(tmp_path):
    """Testing if model notes are written in SBML"""
    path_to_file = join(str(tmp_path), "model_notes.xml")

    # making a minimal cobra model to test notes
    model = cobra.Model("e_coli_core")
    model.notes["Remark"] = "...Model Notes..."
    met = cobra.Metabolite("pyr_c", compartment="c")
    model.add_metabolites([met])
    met.notes["Remark"] = "Note with \n newline"
    rxn = cobra.Reaction("R_ATPM")
    model.add_reactions([rxn])
    rxn.notes["Remark"] = "What about me?"
    model.objective_direction = "max"
    model.objective = rxn
    write_sbml_model(model, path_to_file)

    # reading the model back
    model_after_reading = read_sbml_model(path_to_file)
    met_after_reading = model_after_reading.metabolites.get_by_id("pyr_c")
    reaction_after_reading = model_after_reading.reactions.get_by_id("R_ATPM")

    # checking if notes are written to model
    assert model_after_reading.notes["Remark"] == "...Model Notes..."

    # checking notes for metabolite and reaction
    assert met_after_reading.notes["Remark"] == "Note with \n newline"
    assert reaction_after_reading.notes["Remark"] == "What about me?"
