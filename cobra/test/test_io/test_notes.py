from os.path import join

import pytest

import cobra
from cobra.io import read_sbml_model, write_sbml_model
from cobra.util.solver import set_objective


def test_notes(tmp_path):
    """Testing if model notes are written in SBML"""
    path_to_file = join(str(tmp_path), "model_notes.xml")

    # making a minimal cobra model to test notes
    model = cobra.Model("e_coli_core")
    model.notes["Remark"] = "...Model Notes..."
    met = cobra.Metabolite("pyr_c", compartment="c")
    model.add_metabolites([met])
    met.notes["Remark"] = "I appear."
    rxn = cobra.Reaction("R_ATPM")
    model.add_reactions([rxn])
    rxn.notes["Remark"] = "What about me?"
    model.objective_direction = "max"
    coefficients = {}
    coefficients[rxn] = 1
    set_objective(model, coefficients)
    write_sbml_model(model, path_to_file)

    # reading the model back
    model_after_reading = read_sbml_model(path_to_file)

    # checking if notes are written to model
    assert model_after_reading.notes["Remark"] == "...Model Notes..."
