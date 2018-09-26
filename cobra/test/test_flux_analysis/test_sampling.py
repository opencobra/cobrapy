# -*- coding: utf-8 -*-

"""Test functionalities of flux sampling methods."""

from __future__ import absolute_import

import numpy as np

import pytest
from cobra.core import Metabolite, Model, Reaction
from cobra.flux_analysis.parsimonious import pfba
from cobra.flux_analysis.sampling import ACHRSampler, OptGPSampler, sample


@pytest.fixture(scope="function")
def achr(model):
    """Return ACHRSampler instance for tests."""
    sampler = ACHRSampler(model, thinning=1)
    assert ((sampler.n_warmup > 0) and
            (sampler.n_warmup <= 2 * len(model.variables)))
    assert all(sampler.validate(sampler.warmup) == "v")
    return sampler


@pytest.fixture(scope="function")
def optgp(model):
    """Return OptGPSampler instance for tests."""
    sampler = OptGPSampler(model, processes=1, thinning=1)
    assert ((sampler.n_warmup > 0) and
            (sampler.n_warmup <= 2 * len(model.variables)))
    assert all(sampler.validate(sampler.warmup) == "v")
    return sampler


# Benchmarks
def test_achr_init_benchmark(model, benchmark):
    """Benchmark inital ACHR sampling."""
    benchmark(lambda: ACHRSampler(model))


def test_optgp_init_benchmark(model, benchmark):
    """Benchmark inital OptGP sampling."""
    benchmark(lambda: OptGPSampler(model, processes=2))


def test_achr_sample_benchmark(achr, benchmark):
    """Benchmark ACHR sampling."""
    benchmark(achr.sample, 1)


def test_optgp_sample_benchmark(optgp, benchmark):
    """Benchmark OptGP sampling."""
    benchmark(optgp.sample, 1)


# Tests
def test_single_achr(model):
    """Test ACHR sampling (one sample)."""
    s = sample(model, 10, method="achr")
    assert s.shape == (10, len(model.reactions))


def test_single_optgp(model):
    """Test OptGP sampling (one sample)."""
    s = sample(model, 10, processes=1)
    assert s.shape == (10, len(model.reactions))


def test_multi_optgp(model):
    """Test OptGP sampling (multi sample)."""
    s = sample(model, 10, processes=2)
    assert s.shape == (10, len(model.reactions))


def test_wrong_method(model):
    """Test method intake sanity."""
    with pytest.raises(ValueError):
        sample(model, 1, method="schwupdiwupp")


def test_validate_wrong_sample(achr, model):
    """Test sample correctness."""
    s = achr.sample(10)
    s["hello"] = 1
    with pytest.raises(ValueError):
        achr.validate(s)


def test_fixed_seed(model):
    """Test result of fixed seed for sampling."""
    s1 = sample(model, 1, seed=42)
    s2 = sample(model, 1, seed=42)
    assert np.isclose(s1.TPI[0], s2.TPI[0])


def test_equality_constraint(model):
    """Test equality constraint."""
    model.reactions.ACALD.bounds = (-1.5, -1.5)
    s = sample(model, 10)
    assert np.allclose(s.ACALD, -1.5, atol=1e-6, rtol=0)
    s = sample(model, 10, method="achr")
    assert np.allclose(s.ACALD, -1.5, atol=1e-6, rtol=0)


def test_inequality_constraint(model):
    """Test inequality constraint."""
    co = model.problem.Constraint(
        model.reactions.ACALD.flux_expression, lb=-0.5)
    model.add_cons_vars(co)
    s = sample(model, 10)
    assert all(s.ACALD > -0.5 - 1e-6)
    s = sample(model, 10, method="achr")
    assert all(s.ACALD > -0.5 - 1e-6)


def test_sampling(achr, optgp):
    """Test sampling."""
    s = achr.sample(10)
    assert all(achr.validate(s) == "v")

    s = optgp.sample(10)
    assert all(optgp.validate(s) == "v")


def test_batch_sampling(achr, optgp):
    """Test batch sampling."""
    for b in achr.batch(5, 4):
        assert all(achr.validate(b) == "v")

    for b in optgp.batch(5, 4):
        assert all(optgp.validate(b) == "v")


def test_variables_samples(achr, optgp):
    """Test variable samples."""
    vnames = np.array([v.name for v in achr.model.variables])
    s = achr.sample(10, fluxes=False)
    assert s.shape == (10, achr.warmup.shape[1])
    assert (s.columns == vnames).all()
    assert (achr.validate(s) == "v").all()

    s = optgp.sample(10, fluxes=False)
    assert s.shape == (10, optgp.warmup.shape[1])
    assert (s.columns == vnames).all()
    assert (optgp.validate(s) == "v").all()


def test_inhomogeneous_sanity(model):
    """Test whether inhomogeneous sampling gives approximately the same
    standard deviation as a homogeneous version."""
    model.reactions.ACALD.bounds = (-1.5, -1.5)
    s_inhom = sample(model, 64)

    model.reactions.ACALD.bounds = (-1.5 - 1e-3, -1.5 + 1e-3)
    s_hom = sample(model, 64)

    relative_diff = (s_inhom.std() + 1e-12) / (s_hom.std() + 1e-12)
    assert 0.5 < relative_diff.abs().mean() < 2

    model.reactions.ACALD.bounds = (-1.5, -1.5)
    s_inhom = sample(model, 64, method="achr")

    model.reactions.ACALD.bounds = (-1.5 - 1e-3, -1.5 + 1e-3)
    s_hom = sample(model, 64, method="achr")

    relative_diff = (s_inhom.std() + 1e-12) / (s_hom.std() + 1e-12)
    assert 0.5 < relative_diff.abs().mean() < 2


def test_reproject(optgp):
    """Test reprojection of sampling."""
    s = optgp.sample(10, fluxes=False).values
    proj = np.apply_along_axis(optgp._reproject, 1, s)
    assert all(optgp.validate(proj) == "v")

    s = np.random.rand(10, optgp.warmup.shape[1])
    proj = np.apply_along_axis(optgp._reproject, 1, s)
    assert all(optgp.validate(proj) == "v")


def test_complicated_model():
    """Test a complicated model.

    Difficult model since the online mean calculation is numerically
    unstable so many samples weakly violate the equality constraints.

    """
    model = Model('flux_split')

    reaction1 = Reaction('V1')
    reaction2 = Reaction('V2')
    reaction3 = Reaction('V3')
    reaction1.bounds = (0, 6)
    reaction2.bounds = (0, 8)
    reaction3.bounds = (0, 10)

    A = Metabolite('A')

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


def test_single_point_space(model):
    """Test the reduction of the sampling space to one point."""
    pfba_sol = pfba(model)
    pfba_const = model.problem.Constraint(
        sum(model.variables), ub=pfba_sol.objective_value)
    model.add_cons_vars(pfba_const)
    model.reactions.Biomass_Ecoli_core.lower_bound = \
        pfba_sol.fluxes.Biomass_Ecoli_core
    with pytest.raises(ValueError):
        s = sample(model, 1)
