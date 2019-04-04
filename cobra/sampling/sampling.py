# -*- coding: utf-8 -*-

"""Module implementing flux sampling for cobra models."""

from __future__ import absolute_import, division

import pandas

from cobra.sampling.achr import ACHRSampler
from cobra.sampling.optgp import OptGPSampler


def sample(model, n, method="optgp", thinning=100, processes=1, seed=None):
    """Sample valid flux distributions from a cobra model.

    The function samples valid flux distributions from a cobra model.
    Currently we support two methods:

    1. 'optgp' (default) which uses the OptGPSampler that supports parallel
        sampling [1]_. Requires large numbers of samples to be performant
        (n < 1000). For smaller samples 'achr' might be better suited.

    or

    2. 'achr' which uses artificial centering hit-and-run. This is a single
       process method with good convergence [2]_.

    Parameters
    ----------
    model : cobra.Model
        The model from which to sample flux distributions.
    n : int
        The number of samples to obtain. When using 'optgp' this must be a
        multiple of `processes`, otherwise a larger number of samples will be
        returned.
    method : str, optional
        The sampling algorithm to use.
    thinning : int, optional
        The thinning factor of the generated sampling chain. A thinning of 10
        means samples are returned every 10 steps. Defaults to 100 which in
        benchmarks gives approximately uncorrelated samples. If set to one
        will return all iterates.
    processes : int, optional
        Only used for 'optgp'. The number of processes used to generate
        samples.
    seed : int > 0, optional
        The random number seed to be used. Initialized to current time stamp
        if None.

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
    .. [2] Direction Choice for Accelerated Convergence in Hit-and-Run Sampling
       David E. Kaufman Robert L. Smith
       Operations Research 199846:1 , 84-95

    """

    if method == "optgp":
        sampler = OptGPSampler(model, processes, thinning=thinning, seed=seed)
    elif method == "achr":
        sampler = ACHRSampler(model, thinning=thinning, seed=seed)
    else:
        raise ValueError("method must be 'optgp' or 'achr'!")

    return pandas.DataFrame(columns=[rxn.id for rxn in model.reactions],
                            data=sampler.sample(n))
