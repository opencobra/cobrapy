"""Unit test the ReactionSummary class."""


import pytest

from cobra.summary import ReactionSummary


@pytest.mark.parametrize(
    "reaction_id, expected", [("ACALD", 0.0), ("FUM", 5.06), ("PFK", 7.48),]
)
def test_reaction_summary_flux(model, reaction_id: str, expected: float) -> None:
    """Test that the reported flux in the summary is reasonable."""
    result = ReactionSummary(
        reaction=model.reactions.get_by_id(reaction_id), model=model
    )
    assert result.flux.at[reaction_id, "flux"] == pytest.approx(expected, abs=1e-2)


@pytest.mark.parametrize(
    "reaction_id, min_flux, max_flux",
    [("ACALD", -1.27, 0.0), ("FUM", 0.79, 7.38), ("PFK", 2.58, 16.38)],
)
def test_reaction_summary_flux_fva(
    model, reaction_id: str, min_flux: float, max_flux: float
) -> None:
    """Test that the reported flux ranges in the summary are reasonable."""
    result = ReactionSummary(
        reaction=model.reactions.get_by_id(reaction_id), model=model, fva=0.95
    )
    assert result.flux.at[reaction_id, "minimum"] == pytest.approx(min_flux, abs=1e-2)
    assert result.flux.at[reaction_id, "maximum"] == pytest.approx(max_flux, abs=1e-2)
