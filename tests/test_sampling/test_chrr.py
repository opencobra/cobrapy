"""Test functionalities of flux sampling methods."""

import pytest

from cobra.core import Metabolite, Model, Reaction
from cobra.sampling import hopsy_is_available


@pytest.mark.skipif(not hopsy_is_available, reason="hopsy not installed")
def test_complicated_model_chrr() -> None:
    """Test a complicated model with CHRR.

    Difficult model since the online mean calculation is numerically
    unstable, so many samples weakly violate the equality constraints.

    """
    from cobra.sampling import HopsySampler

    model = Model("flux_split")

    reaction1 = Reaction("V1")
    reaction2 = Reaction("V2")
    reaction3 = Reaction("V3")
    reaction1.bounds = (0, 6)
    reaction2.bounds = (0, 8)
    reaction3.bounds = (0, 10)

    A = Metabolite("A")

    reaction1.add_metabolites({A: -1})
    reaction2.add_metabolites({A: -1})
    reaction3.add_metabolites({A: 1})

    model.add_reactions([reaction1, reaction2, reaction3])

    chrr = HopsySampler(model, seed=42)
    chrr_samples = chrr.sample(100)

    assert any(chrr_samples.corr().abs() < 1.0)
    assert sum(chrr.validate(chrr_samples) == "v") > 95
