"""Unit test the model summary."""


import pytest

from cobra.flux_analysis import flux_variability_analysis, pfba
from cobra.summary import ModelSummary


def test_model_summary_interface(model, opt_solver):
    """Test that a summary can be created successfully."""
    model.solver = opt_solver
    ModelSummary(
        model=model,
    )
    ModelSummary(
        model=model,
        solution=pfba(model),
    )
    ModelSummary(
        model=model,
        fva=0.95,
    )
    ModelSummary(
        model=model,
        fva=flux_variability_analysis(
            model, reaction_list=["CYTBD", "NADH16", "SUCDi"]
        ),
    )


def test_model_summary_to_frame(model, opt_solver):
    """Test that the summary's method ``to_frame`` can be called."""
    model.solver = opt_solver
    summary = model.summary()
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
def test_model_summary_to_string(model, opt_solver, kwargs):
    """Test that the summary's method ``to_string`` can be called."""
    model.solver = opt_solver
    summary = model.summary()
    summary.to_string(**kwargs)


@pytest.mark.parametrize(
    "kwargs", [{}, {"names": True}, {"float_format": ".1f"}, {"threshold": 1e-2}]
)
def test_model_summary_to_html(model, opt_solver, kwargs):
    """Test that the summary's method ``to_html`` can be called."""
    model.solver = opt_solver
    summary = model.summary()
    summary.to_html(**kwargs)


def test_model_summary_to_frame_previous_solution(model, opt_solver):
    """Test that the summary correctly uses an existing solution."""
    model.solver = opt_solver
    solution = pfba(model)
    rxn = model.reactions.EX_glc__D_e
    solution.fluxes[rxn.id] = -123
    summary = model.summary(solution=solution)
    assert summary.to_frame().at[rxn.id, "flux"] == 123


def test_model_summary_flux(model, opt_solver):
    """Test that the summary has expected fluxes."""
    model.solver = opt_solver
    summary = model.summary()
    assert summary.uptake_flux.at["EX_o2_e", "flux"] == pytest.approx(21.8, abs=1e-02)
    assert summary.secretion_flux.at["EX_h2o_e", "flux"] == pytest.approx(
        -29.18, abs=1e-02
    )


def test_model_summary_fva(model, opt_solver):
    """Test that the exchange summary is within expected bounds."""
    model.solver = opt_solver
    summary = model.summary(fva=0.95)
    assert summary.to_frame().at["EX_o2_e", "flux"] == pytest.approx(21.8, abs=1e-02)
    assert summary.to_frame().at["EX_o2_e", "minimum"] == pytest.approx(19.9, abs=1e-02)
    assert summary.to_frame().at["EX_o2_e", "maximum"] == pytest.approx(
        23.71, abs=1e-02
    )

    assert summary.to_frame().at["EX_h2o_e", "flux"] == pytest.approx(-29.18, abs=1e-02)
    assert summary.to_frame().at["EX_h2o_e", "minimum"] == pytest.approx(
        -30.72, abs=1e-02
    )
    assert summary.to_frame().at["EX_h2o_e", "maximum"] == pytest.approx(-25, abs=1e-02)


def test_model_summary_flux_in_context(model, opt_solver):
    """Test that the model summary inside and outside of a context are equal."""
    model.solver = opt_solver
    copy = model.copy()
    with model:
        model.reactions.EX_glc__D_e.bounds = (0, 1000)
        model.reactions.EX_ac_e.bounds = (-10, 1000)
        context_summary = model.summary()

    copy.reactions.EX_glc__D_e.bounds = (0, 1000)
    copy.reactions.EX_ac_e.bounds = (-10, 1000)
    outside_summary = copy.summary()

    assert context_summary.to_frame()["flux"].values == pytest.approx(
        outside_summary.to_frame()["flux"].values, abs=model.tolerance
    )
