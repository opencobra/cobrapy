# -*- coding: utf-8 -*-

"""Test functionalities of ACHRSampler."""

from __future__ import absolute_import

import numpy as np
import pytest

from cobra.sampling import ACHRSampler


def test_achr_init_benchmark(model, benchmark):
    """Benchmark inital ACHR sampling."""

    benchmark(lambda: ACHRSampler(model))


def test_achr_sample_benchmark(achr, benchmark):
    """Benchmark ACHR sampling."""

    benchmark(achr.sample, 1)


def test_validate_wrong_sample(achr, model):
    """Test sample correctness."""

    s = achr.sample(10)
    s["hello"] = 1

    with pytest.raises(ValueError):
        achr.validate(s)


def test_sampling(achr):
    """Test sampling."""

    s = achr.sample(10)
    assert all(achr.validate(s) == "v")


def test_batch_sampling(achr):
    """Test batch sampling."""

    for b in achr.batch(5, 4):
        assert all(achr.validate(b) == "v")


def test_variables_samples(achr):
    """Test variable samples."""

    vnames = np.array([v.name for v in achr.model.variables])
    s = achr.sample(10, fluxes=False)
    assert s.shape == (10, achr.warmup.shape[1])
    assert (s.columns == vnames).all()
    assert (achr.validate(s) == "v").all()
