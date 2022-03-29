"""Test functionalities of ACHRSampler."""

from typing import TYPE_CHECKING, Callable

import numpy as np
import pytest

from cobra.sampling import ACHRSampler


if TYPE_CHECKING:
    from cobra import Model


def test_achr_init_benchmark(model: "Model", benchmark: Callable) -> None:
    """Benchmark inital ACHR sampling."""
    benchmark(lambda: ACHRSampler(model))


def test_achr_sample_benchmark(achr: ACHRSampler, benchmark: Callable) -> None:
    """Benchmark ACHR sampling."""
    benchmark(achr.sample, 1)


def test_validate_wrong_sample(achr: ACHRSampler, model: "Model") -> None:
    """Test sample correctness."""
    s = achr.sample(10)
    s["hello"] = 1

    with pytest.raises(ValueError):
        achr.validate(s)


def test_sampling(achr: ACHRSampler) -> None:
    """Test sampling."""
    s = achr.sample(10)
    assert all(achr.validate(s) == "v")


def test_batch_sampling(achr: ACHRSampler) -> None:
    """Test batch sampling."""
    for b in achr.batch(5, 4):
        assert all(achr.validate(b) == "v")


def test_variables_samples(achr: ACHRSampler) -> None:
    """Test variable samples."""
    vnames = np.array([v.name for v in achr.model.variables])
    s = achr.sample(10, fluxes=False)
    assert s.shape == (10, achr.warmup.shape[1])
    assert (s.columns == vnames).all()
    assert (achr.validate(s) == "v").all()
