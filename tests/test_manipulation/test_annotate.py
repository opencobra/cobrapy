"""Test functionalities of model component annotation functions."""

from cobra.core import Model, Reaction
from cobra.manipulation import add_SBO


def test_sbo_annotation(model: Model) -> None:
    """Test SBO annotation function."""
    rxns = model.reactions
    rxns.EX_o2_e.annotation.clear()
    fake_DM = Reaction("DM_h_c")
    model.add_reactions([fake_DM])
    fake_DM.add_metabolites({model.metabolites.get_by_id("h_c"): -1})
    # this exchange will be set wrong. The function should not overwrite
    # an existing SBO annotation
    rxns.get_by_id("EX_h_e").annotation["sbo"] = "SBO:0000628"
    add_SBO(model)
    assert rxns.EX_o2_e.annotation["sbo"] == "SBO:0000627"
    assert rxns.DM_h_c.annotation["sbo"] == "SBO:0000628"
    assert rxns.EX_h_e.annotation["sbo"] == "SBO:0000628"
