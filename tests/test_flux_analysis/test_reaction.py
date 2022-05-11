"""Test assessing functions in flux_analysis.reaction."""

from typing import List

from cobra.core import Metabolite, Model, Reaction
from cobra.flux_analysis.reaction import assess


def test_assess(model: Model, all_solvers: List[str]) -> None:
    """Test assess functions."""
    with model:
        assert assess(model, model.reactions.GLCpts, solver=all_solvers) is True
        pyr = model.metabolites.pyr_c
        a = Metabolite("a")
        b = Metabolite("b")
        model.add_metabolites([a, b])
        pyr_a2b = Reaction("pyr_a2b")
        pyr_a2b.add_metabolites({pyr: -1, a: -1, b: 1})
        model.add_reactions([pyr_a2b])
        res = assess(model, pyr_a2b, 0.01, solver=all_solvers)
        expected = {
            "precursors": {a: {"required": 0.01, "produced": 0.0}},
            "products": {b: {"required": 0.01, "capacity": 0.0}},
        }
        assert res == expected
