"""Unit test the MetaboliteSummary class."""


import pytest

from cobra.flux_analysis import flux_variability_analysis, pfba
from cobra.summary import MetaboliteSummary


def test_metabolite_summary_interface(model, opt_solver):
    """Test that a summary can be created successfully."""
    model.solver = opt_solver
    metabolite = model.metabolites.get_by_id("q8_c")
    MetaboliteSummary(
        metabolite=metabolite,
        model=model,
    )
    MetaboliteSummary(
        metabolite=metabolite,
        model=model,
        solution=pfba(model),
    )
    MetaboliteSummary(
        metabolite=metabolite,
        model=model,
        fva=0.95,
    )
    MetaboliteSummary(
        metabolite=metabolite,
        model=model,
        fva=flux_variability_analysis(
            model, reaction_list=["CYTBD", "NADH16", "SUCDi"]
        ),
    )


def test_metabolite_summary_to_frame(model, opt_solver):
    """Test that the summary's method ``to_frame`` can be called."""
    model.solver = opt_solver
    summary = model.metabolites.get_by_id("atp_c").summary()
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
def test_metabolite_summary_to_string(model, opt_solver, kwargs):
    """Test that the summary's method ``to_string`` can be called."""
    model.solver = opt_solver
    summary = model.metabolites.get_by_id("atp_c").summary()
    summary.to_string(**kwargs)


@pytest.mark.parametrize(
    "kwargs", [{}, {"names": True}, {"float_format": ".1f"}, {"threshold": 1e-2}]
)
def test_metabolite_summary_to_html(model, opt_solver, kwargs):
    """Test that the summary's method ``to_html`` can be called."""
    model.solver = opt_solver
    summary = model.metabolites.get_by_id("atp_c").summary()
    summary.to_html(**kwargs)


def test_q8_producing_summary(model, opt_solver):
    """Test that the production summary of q8 is accurate."""
    model.solver = opt_solver
    summary = model.metabolites.get_by_id("q8_c").summary()

    assert summary.producing_flux.at["CYTBD", "percent"] == 1
    assert summary.producing_flux.at["CYTBD", "flux"] == pytest.approx(43.6, abs=1e-2)


def test_q8_consuming_summary(model, opt_solver):
    """Test that the consumption summary of q8 is accurate."""
    model.solver = opt_solver
    summary = model.metabolites.get_by_id("q8_c").summary()

    assert summary.consuming_flux.at["NADH16", "percent"] == pytest.approx(
        0.8838, abs=1e-4
    )
    assert summary.consuming_flux.at["NADH16", "flux"] == pytest.approx(
        -38.54, abs=1e-2
    )

    assert summary.consuming_flux.at["SUCDi", "percent"] == pytest.approx(
        0.1162, abs=1e-4
    )
    assert summary.consuming_flux.at["SUCDi", "flux"] == pytest.approx(-5.06, abs=1e-2)


def test_fdp_production_with_fva(model, opt_solver):
    """Test that the production summary of fdp is within expected bounds."""
    model.solver = opt_solver
    summary = model.metabolites.get_by_id("fdp_c").summary(fva=0.99)
    assert summary.producing_flux.at["PFK", "percent"] == 1
    assert summary.producing_flux.at["PFK", "flux"] == pytest.approx(7.48, abs=1e-2)
    assert summary.producing_flux.at["PFK", "minimum"] == pytest.approx(6.17, abs=1e-2)
    assert summary.producing_flux.at["PFK", "maximum"] == pytest.approx(9.26, abs=1e-2)


@pytest.mark.parametrize("metabolite_id", ["q8_c", "fdp_c", "atp_c"])
def test_metabolite_summary_flux_in_context(model, opt_solver, metabolite_id: str):
    """Test that the metabolite summary inside and outside of a context are equal."""
    model.solver = opt_solver
    with model:
        context_summary = model.metabolites.get_by_id(metabolite_id).summary()
    outside_summary = model.metabolites.get_by_id(metabolite_id).summary()

    assert context_summary.to_frame()["flux"].values == pytest.approx(
        outside_summary.to_frame()["flux"].values, abs=model.tolerance
    )
