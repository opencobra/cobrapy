"""Provide a wrapper function for performing flux sampling of cobra models."""

import logging
from typing import TYPE_CHECKING, Optional

import pandas as pd

from .achr import ACHRSampler
from .hopsy import hopsy_is_available
from .optgp import OptGPSampler


logger = logging.getLogger(__name__)


if hopsy_is_available:
    from .hopsy import HopsySampler


if TYPE_CHECKING:
    from cobra import Model


def sample(
    model: "Model",
    n: int,
    method: str = "auto",
    thinning: int = 100,
    processes: int = 1,
    seed: Optional[int] = None,
) -> pd.DataFrame:
    """Sample valid flux distributions from a cobra model.

    Currently, three methods are supported:

    1. 'chrr' which uses a the coordinate Hit-and-Run algorithm
        with rounding ([1]_), which has been shown to often outperform OptGP
        and ACHR, refer [2]_, [3]_. The rounding transformation is performed
        in an offline manner which inflicts a certain base cost, however
        sampling is performed using the C++-based library ``hopsy`` ([4]_)
        and is, thus, blazingly fast. Moreover, this method supports
        parallel sampling.

    2. 'optgp'  which uses the OptGPSampler that supports parallel
        sampling. Requires large numbers of samples to be performant
        (`n` > 1000). For smaller samples, 'achr' might be better suited.
        For details, refer [5]_ .

    3. 'achr' which uses artificial centering hit-and-run. This is a single
       process method with good convergence. For details, refer [6]_ .

    A fourth 'auto' option is available, which tries to choose 'chrr' if ``hopsy``
    is available. This is the default option.

    Parameters
    ----------
    model : cobra.Model
        The model from which to sample flux distributions.
    n : int
        The number of samples to obtain. When using 'optgp', this must be a
        multiple of `processes`, otherwise a larger number of samples will
        be returned.
    method : {"auto", "chrr", "optgp", "achr"}, optional
        The sampling algorithm to use (default "auto").
    thinning : int, optional
        The thinning factor of the generated sampling chain. A thinning of
        10 means samples are returned every 10 steps. Defaults to 100 which
        in benchmarks gives approximately uncorrelated samples. If set to 1
        will return all iterates (default 100).
    processes : int, optional
        Only used for 'optgp' and 'chrr'. The number of processes
        used to generate samples (default 1).
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
    .. [1] Hulda S Haraldsdóttir, Ben Cousins, Ines Thiele, Ronan M.T Fleming,
       Santosh Vempala,
       CHRR: coordinate hit-and-run with rounding for uniform sampling of
       constraint-based models,
       Bioinformatics, Volume 33, Issue 11, June 2017, Pages 1741–1743,
       https://doi.org/10.1093/bioinformatics/btx052

    .. [2] Herrmann, H.A., Dyson, B.C., Vass, L. et al.
       Flux sampling is a powerful tool to study metabolism under changing
       environmental conditions.
       npj Syst Biol Appl 5, 32 (2019).
       https://doi.org/10.1038/s41540-019-0109-0

    .. [3] Fallahi S, Skaug HJ, Alendal G (2020)
       A comparison of Monte Carlo sampling methods for
       metabolic network models.
       PLoS ONE 15(7): e0235393.
       https://doi.org/10.1371/journal.pone.0235393

    .. [4] Richard D Paul, Johann F Jadebeck, Anton Stratmann, Wolfgang Wiechert,
       Katharina Nöh,
       hopsy — a methods marketplace for convex polytope sampling in Python,
       Bioinformatics, Volume 40, Issue 7, July 2024, btae430,
       https://doi.org/10.1093/bioinformatics/btae430

    .. [5] Megchelenbrink W, Huynen M, Marchiori E (2014)
       optGpSampler: An Improved Tool for Uniformly Sampling the Solution-Space
       of Genome-Scale Metabolic Networks.
       PLoS ONE 9(2): e86587.
       https://doi.org/10.1371/journal.pone.0086587

    .. [6] Direction Choice for Accelerated Convergence in Hit-and-Run Sampling
       David E. Kaufman, Robert L. Smith
       Operations Research 199846:1 , 84-95
       https://doi.org/10.1287/opre.46.1.84

    """
    if method == "auto":
        if hopsy_is_available:
            method = "chrr"
        else:
            method = "optgp"
            logger.warn(
                "hopsy is not available in your environment."
                "Falling back to 'optgp' sampler.",
                stacklevel=2,
            )

    if method == "optgp":
        sampler = OptGPSampler(model, processes=processes, thinning=thinning, seed=seed)
    elif method == "achr":
        sampler = ACHRSampler(model, thinning=thinning, seed=seed)
    elif method == "chrr":
        sampler = HopsySampler(
            model, processes=processes, thinning=thinning, seed=seed, rounding=True
        )
    else:
        raise ValueError(
            f'Invalid value: "{method}" for method used. '
            'The value must be "optgp", "achr" or "chrr".'
        )

    return sampler.sample(n)
