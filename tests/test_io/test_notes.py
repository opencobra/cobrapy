"""Test proper reading of SBML notes."""

from pathlib import Path

from cobra import Metabolite, Model, Reaction
from cobra.io import read_sbml_model, write_sbml_model


def test_notes(tmp_path: Path) -> None:
    """Test if model notes are written in SBML.

    Parameters
    ----------
    tmp_path : pathlib.Path
        The path to the temporary test assets store.

    """
    path_to_file = tmp_path / "model_notes.xml"

    # making a minimal cobra model to test notes
    model = Model("test_notes_model")
    model.notes["Remark"] = "...Model Notes..."
    met = Metabolite("pyr_c", compartment="c")
    model.add_metabolites([met])
    met.notes["Remark"] = "Note with \n newline"
    rxn = Reaction("R_ATPM")
    model.add_reactions([rxn])
    rxn.notes["Remark"] = "What about me?"
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
