"""Define module level fixtures."""

from typing import TYPE_CHECKING

import pytest

from cobra.sampling import ACHRSampler


if TYPE_CHECKING:
    from cobra import Model


@pytest.fixture(scope="function")
def achr(model: "Model") -> ACHRSampler:
    """Return ACHRSampler instance for tests."""
    sampler = ACHRSampler(model, thinning=1)
    assert (sampler.n_warmup > 0) and (sampler.n_warmup <= 2 * len(model.variables))
    assert all(sampler.validate(sampler.warmup) == "v")

    return sampler
