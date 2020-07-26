# -*- coding: utf-8 -*-

"""Define module level fixtures."""

from __future__ import absolute_import

import pytest

from cobra.sampling import ACHRSampler


@pytest.fixture(scope="function")
def achr(model):
    """Return ACHRSampler instance for tests."""

    sampler = ACHRSampler(model, thinning=1)
    assert (sampler.n_warmup > 0) and (sampler.n_warmup <= 2 * len(model.variables))
    assert all(sampler.validate(sampler.warmup) == "v")

    return sampler
