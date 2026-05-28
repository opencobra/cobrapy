"""Test functionalities of flux sampling methods."""

import os

import numpy as np
import pytest

from cobra.core import Metabolite, Model, Reaction
from cobra.flux_analysis.parsimonious import pfba
from cobra.sampling import ACHRSampler, OptGPSampler, hopsy_is_available, sample


def test_single_achr(model: Model) -> None:
    """Test ACHR sampling (one sample)."""
    s = sample(model, 10, method="achr")
    assert s.shape == (10, len(model.reactions))


@pytest.mark.skipif(not hopsy_is_available, reason="hopsy not installed")
def test_single_chrr(model: Model) -> None:
    """Test CHRR sampling (one sample)."""
    s = sample(model, 10, processes=1, method="chrr")
    assert s.shape == (10, len(model.reactions))


def test_single_optgp(model: Model) -> None:
    """Test OptGP sampling (one sample)."""
    s = sample(model, 10, processes=1, method="optgp")
    assert s.shape == (10, len(model.reactions))


@pytest.mark.skipif(
    "SKIP_MP" in os.environ or not hopsy_is_available,
    reason="unsafe for parallel execution or hopsy not installed",
)
def test_multi_chrr(model: Model) -> None:  # pragma: no cover
    """Test CHRR sampling (multi sample)."""
    s = sample(model, 10, processes=2, method="chrr")
    assert s.shape == (10, len(model.reactions))


@pytest.mark.skipif("SKIP_MP" in os.environ, reason="unsafe for parallel execution")
def test_multi_optgp(model: Model) -> None:  # pragma: no cover
    """Test OptGP sampling (multi sample)."""
    s = sample(model, 10, processes=2, method="optgp")
    assert s.shape == (10, len(model.reactions))


def test_wrong_method(model: Model) -> None:
    """Test method intake sanity."""
    with pytest.raises(ValueError):
        sample(model, 1, method="schwupdiwupp")


def test_fixed_seed(model: Model) -> None:
    """Test result of fixed seed for sampling."""
    s1 = sample(model, 1, seed=42, method="achr")
    s2 = sample(model, 1, seed=42, method="achr")
    assert np.isclose(s1.TPI[0], s2.TPI[0])

    s1 = sample(model, 1, seed=42, method="optgp")
    s2 = sample(model, 1, seed=42, method="optgp")
    assert np.isclose(s1.TPI[0], s2.TPI[0])


@pytest.mark.skipif(not hopsy_is_available, reason="hopsy not installed")
def test_fixed_seed_chrr(model: Model) -> None:
    """Test result of fixed seed for CHRR sampling."""
    s1 = sample(model, 1, seed=42, method="chrr")
    s2 = sample(model, 1, seed=42, method="chrr")
    assert np.isclose(s1.TPI[0], s2.TPI[0])


def test_equality_constraint(model: Model) -> None:
    """Test equality constraint."""
    model.reactions.ACALD.bounds = (-1.5, -1.5)

    s = sample(model, 10, method="achr")
    assert np.allclose(s.ACALD, -1.5, atol=1e-6, rtol=0)

    s = sample(model, 10, method="optgp")
    assert np.allclose(s.ACALD, -1.5, atol=1e-6, rtol=0)


@pytest.mark.skipif(not hopsy_is_available, reason="hopsy not installed")
def test_equality_constraint_chrr(model: Model) -> None:
    """Test equality constraint with CHRR."""
    model.reactions.ACALD.bounds = (-1.5, -1.5)

    s = sample(model, 10, method="chrr")
    assert np.allclose(s.ACALD, -1.5, atol=1e-6, rtol=0)


def test_inequality_constraint(model: Model) -> None:
    """Test inequality constraint."""
    co = model.problem.Constraint(model.reactions.ACALD.flux_expression, lb=-0.5)
    model.add_cons_vars(co)

    s = sample(model, 10, method="achr")
    assert all(s.ACALD > -0.5 - 1e-6)

    s = sample(model, 10, method="optgp")
    assert all(s.ACALD > -0.5 - 1e-6)


@pytest.mark.skipif(not hopsy_is_available, reason="hopsy not installed")
def test_inequality_constraint_chrr(model: Model) -> None:
    """Test inequality constraint with CHRR."""
    co = model.problem.Constraint(model.reactions.ACALD.flux_expression, lb=-0.5)
    model.add_cons_vars(co)

    s = sample(model, 10, method="chrr")
    assert all(s.ACALD > -0.5 - 1e-6)


def test_inhomogeneous_sanity(model: Model) -> None:
    """Test standard deviation between inhomogeneous and homogeneous sampling."""
    model.reactions.ACALD.bounds = (-1.5, -1.5)
    s_inhom = sample(model, 64, method="achr")

    model.reactions.ACALD.bounds = (-1.5 - 1e-3, -1.5 + 1e-3)
    s_hom = sample(model, 64, method="achr")

    relative_diff = (s_inhom.std() + 1e-12) / (s_hom.std() + 1e-12)
    assert 0.5 < relative_diff.abs().mean() < 2

    model.reactions.ACALD.bounds = (-1.5, -1.5)
    s_inhom = sample(model, 64, method="optgp")

    model.reactions.ACALD.bounds = (-1.5 - 1e-3, -1.5 + 1e-3)
    s_hom = sample(model, 64, method="optgp")

    relative_diff = (s_inhom.std() + 1e-12) / (s_hom.std() + 1e-12)
    assert 0.5 < relative_diff.abs().mean() < 2


@pytest.mark.skipif(not hopsy_is_available, reason="hopsy not installed")
def test_inhomogeneous_sanity_chrr(model: Model) -> None:
    """Test standard deviation between inhomogeneous and homogeneous CHRR sampling."""
    model.reactions.ACALD.bounds = (-1.5, -1.5)
    s_inhom = sample(model, 64, method="chrr")

    model.reactions.ACALD.bounds = (-1.5 - 1e-3, -1.5 + 1e-3)
    s_hom = sample(model, 64, method="chrr")

    relative_diff = (s_inhom.std() + 1e-12) / (s_hom.std() + 1e-12)
    assert 0.5 < relative_diff.abs().mean() < 2


def test_complicated_model() -> None:
    """Test a complicated model.

    Difficult model since the online mean calculation is numerically
    unstable, so many samples weakly violate the equality constraints.

    """
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

    optgp = OptGPSampler(model, 1, seed=42)
    achr = ACHRSampler(model, seed=42)

    optgp_samples = optgp.sample(100)
    achr_samples = achr.sample(100)

    assert any(optgp_samples.corr().abs() < 1.0)
    assert any(achr_samples.corr().abs() < 1.0)

    # > 95% are valid
    assert sum(optgp.validate(optgp_samples) == "v") > 95
    assert sum(achr.validate(achr_samples) == "v") > 95


def test_single_point_space(model: Model) -> None:
    """Test the reduction of the sampling space to one point."""
    pfba_sol = pfba(model)
    pfba_const = model.problem.Constraint(
        sum(model.variables), ub=pfba_sol.objective_value
    )
    model.add_cons_vars(pfba_const)
    model.reactions.Biomass_Ecoli_core.lower_bound = pfba_sol.fluxes.Biomass_Ecoli_core

    with pytest.raises(ValueError):
        sample(model, 1)
