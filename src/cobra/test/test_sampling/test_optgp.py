"""Test functionalities of OptGPSampler."""

import os
from typing import TYPE_CHECKING, Callable

import numpy as np
import pytest

from cobra.sampling import OptGPSampler


if TYPE_CHECKING:
    from cobra import Model
    from cobra.sampling import ACHRSampler


@pytest.fixture(scope="function")
def optgp(model: "Model") -> OptGPSampler:
    """Return OptGPSampler instance for tests."""
    sampler = OptGPSampler(model, processes=1, thinning=1)
    assert (sampler.n_warmup > 0) and (sampler.n_warmup <= 2 * len(model.variables))
    assert all(sampler.validate(sampler.warmup) == "v")

    return sampler


@pytest.mark.skipif("SKIP_MP" in os.environ, reason="unsafe for parallel execution")
def test_optgp_init_benchmark(model: "Model", benchmark: Callable) -> None:
    """Benchmark inital OptGP sampling."""
    benchmark(lambda: OptGPSampler(model, processes=2))


def test_optgp_sample_benchmark(optgp: "Model", benchmark: Callable) -> None:
    """Benchmark OptGP sampling."""
    benchmark(optgp.sample, 1)


def test_sampling(optgp: OptGPSampler) -> None:
    """Test sampling."""
    s = optgp.sample(10)
    assert all(optgp.validate(s) == "v")


def test_batch_sampling(optgp: OptGPSampler) -> None:
    """Test batch sampling."""
    for b in optgp.batch(5, 4):
        assert all(optgp.validate(b) == "v")


def test_variables_samples(achr: "ACHRSampler", optgp: OptGPSampler) -> None:
    """Test variable samples."""
    vnames = np.array([v.name for v in achr.model.variables])
    s = optgp.sample(10, fluxes=False)
    assert s.shape == (10, optgp.warmup.shape[1])
    assert (s.columns == vnames).all()
    assert (optgp.validate(s) == "v").all()


def test_reproject(optgp: OptGPSampler) -> None:
    """Test reprojection of sampling."""
    s = optgp.sample(10, fluxes=False).values
    proj = np.apply_along_axis(optgp._reproject, 1, s)
    assert all(optgp.validate(proj) == "v")

    s = np.random.rand(10, optgp.warmup.shape[1])
    proj = np.apply_along_axis(optgp._reproject, 1, s)
    assert all(optgp.validate(proj) == "v")
