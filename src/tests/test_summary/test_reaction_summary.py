"""Unit test the ReactionSummary class."""


import pytest

from cobra.flux_analysis import flux_variability_analysis, pfba
from cobra.summary import ReactionSummary


def test_reaction_summary_interface(model, opt_solver):
    """Test that a summary can be created successfully."""
    model.solver = opt_solver
    reaction = model.reactions.get_by_id("FUM")
    ReactionSummary(
        reaction=reaction,
        model=model,
    )
    ReactionSummary(
        reaction=reaction,
        model=model,
        solution=pfba(model),
    )
    ReactionSummary(
        reaction=reaction,
        model=model,
        fva=0.95,
    )
    ReactionSummary(
        reaction=reaction,
        model=model,
        fva=flux_variability_analysis(model, reaction_list=["FUM"]),
    )


def test_reaction_summary_to_frame(model, opt_solver):
    """Test that the summary's method ``to_frame`` can be called."""
    model.solver = opt_solver
    summary = model.reactions.get_by_id("FUM").summary()
    summary.to_frame()


@pytest.mark.parametrize(
    "kwargs",
    [
        {},
        {"names": True},
        {"float_format": ".1f"},
        {"threshold": 1e-2},
        {"column_width": 20},
    ],
)
def test_reaction_summary_to_string(model, opt_solver, kwargs):
    """Test that the summary's method ``to_string`` can be called."""
    model.solver = opt_solver
    summary = model.reactions.get_by_id("FUM").summary()
    summary.to_string(**kwargs)


@pytest.mark.parametrize(
    "kwargs", [{}, {"names": True}, {"float_format": ".1f"}, {"threshold": 1e-2}]
)
def test_reaction_summary_to_html(model, opt_solver, kwargs):
    """Test that the summary's method ``to_html`` can be called."""
    model.solver = opt_solver
    summary = model.reactions.get_by_id("FUM").summary()
    summary.to_html(**kwargs)


@pytest.mark.parametrize(
    "reaction_id, expected", [("ACALD", 0.0), ("FUM", 5.06), ("PFK", 7.48)]
)
def test_reaction_summary_flux(model, reaction_id: str, expected: float) -> None:
    """Test that the reported flux in the summary is reasonable."""
    result = ReactionSummary(
        reaction=model.reactions.get_by_id(reaction_id), model=model
    )
    assert result.to_frame().at[reaction_id, "flux"] == pytest.approx(
        expected, abs=1e-2
    )


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
    assert result.to_frame().at[reaction_id, "minimum"] == pytest.approx(
        min_flux, abs=1e-2
    )
    assert result.to_frame().at[reaction_id, "maximum"] == pytest.approx(
        max_flux, abs=1e-2
    )


@pytest.mark.parametrize("reaction_id", ["ACALD", "FUM", "PFK"])
def test_reaction_summary_flux_in_context(model, reaction_id: str) -> None:
    """Test that the reaction summary inside and outside of a context are equal."""
    with model:
        context_summary = model.reactions.get_by_id(reaction_id).summary()
    outside_summary = model.reactions.get_by_id(reaction_id).summary()

    assert context_summary.to_frame()["flux"].values == pytest.approx(
        outside_summary.to_frame()["flux"].values, abs=model.tolerance
    )
