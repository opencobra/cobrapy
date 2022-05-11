"""Provide a wrapper function for performing flux sampling of cobra models."""

from typing import TYPE_CHECKING, Optional

import pandas as pd

from .achr import ACHRSampler
from .optgp import OptGPSampler


if TYPE_CHECKING:
    from cobra import Model


def sample(
    model: "Model",
    n: int,
    method: str = "optgp",
    thinning: int = 100,
    processes: int = 1,
    seed: Optional[int] = None,
) -> pd.DataFrame:
    """Sample valid flux distributions from a cobra model.

    Currently, two methods are supported:

    1. 'optgp' (default) which uses the OptGPSampler that supports parallel
        sampling. Requires large numbers of samples to be performant
        (`n` > 1000). For smaller samples, 'achr' might be better suited.
        For details, refer [1]_ .

    2. 'achr' which uses artificial centering hit-and-run. This is a single
       process method with good convergence. For details, refer [2]_ .

    Parameters
    ----------
    model : cobra.Model
        The model from which to sample flux distributions.
    n : int
        The number of samples to obtain. When using 'optgp', this must be a
        multiple of `processes`, otherwise a larger number of samples will
        be returned.
    method : {"optgp", "achr"}, optional
        The sampling algorithm to use (default "optgp").
    thinning : int, optional
        The thinning factor of the generated sampling chain. A thinning of
        10 means samples are returned every 10 steps. Defaults to 100 which
        in benchmarks gives approximately uncorrelated samples. If set to 1
        will return all iterates (default 100).
    processes : int, optional
        Only used for 'optgp'. The number of processes used to generate
        samples (default 1).
    seed : int > 0, optional
        Sets the random number seed. Initialized to the current time stamp
        if None (default None).

    Returns
    -------
    pandas.DataFrame
        The generated flux samples. Each row corresponds to a sample of the
        fluxes and the columns are the reactions.

    Notes
    -----
    The samplers have a correction method to ensure equality feasibility for
    long-running chains, however this will only work for homogeneous models,
    meaning models with no non-zero fixed variables or constraints (
    right-hand side of the equalities are zero).

    References
    ----------
    .. [1] Megchelenbrink W, Huynen M, Marchiori E (2014)
       optGpSampler: An Improved Tool for Uniformly Sampling the Solution-Space
       of Genome-Scale Metabolic Networks.
       PLoS ONE 9(2): e86587.
       https://doi.org/10.1371/journal.pone.0086587

    .. [2] Direction Choice for Accelerated Convergence in Hit-and-Run Sampling
       David E. Kaufman, Robert L. Smith
       Operations Research 199846:1 , 84-95
       https://doi.org/10.1287/opre.46.1.84

    """
    if method == "optgp":
        sampler = OptGPSampler(model, processes=processes, thinning=thinning, seed=seed)
    elif method == "achr":
        sampler = ACHRSampler(model, thinning=thinning, seed=seed)
    else:
        raise ValueError(
            f'Invalid value: "{method}" for method used. '
            'The value must be "optgp" or "achr".'
        )

    return pd.DataFrame(
        columns=[rxn.id for rxn in model.reactions], data=sampler.sample(n)
    )
