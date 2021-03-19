"""Provide low-level sampling stepper functions and helpers."""

import logging
from typing import Optional

import numpy as np

from .hr_sampler import HRSampler


logger = logging.getLogger(__name__)


# Maximum number of retries for sampling
MAX_TRIES = 100


def step(
    sampler: HRSampler,
    x: np.ndarray,
    delta: np.ndarray,
    fraction: Optional[float] = None,
    tries: int = 0,
) -> np.ndarray:
    """Sample a new feasible point from the point `x` in direction `delta`.

    This is the low-level sampling stepper for samplers derived
    from `HRSampler`. Currently, it's used by `ACHRSampler` and
    `OptGPSampler`.

    It's declared outside of the base sampling class to facilitate use of
    multiprocessing.

    Parameters
    ----------
    sampler : cobra.sampling.HRSampler
        The sampler to sample a step for.
    x : np.array
        A point in the sampling region.
    delta : np.array
        The direction to take the step in.
    fraction : float, optional
        A float controlling the part of alpha difference to contribute to
        the fraction of `delta` (default None). If None, alpha is obtained
        from a normal distribution.
    tries : int, optional
        Total number of tries (default 0).

    Returns
    -------
    np.array
        The new numpy array obtained after a step of sampling.

    Raises
    ------
    RuntimeError
        If `tries` exceeds `MAX_TRIES`.

    """
    prob = sampler.problem
    valid = (np.abs(delta) > sampler.feasibility_tol) & np.logical_not(
        prob.variable_fixed
    )

    # permissible alphas for staying in variable bounds
    valphas = ((1.0 - sampler.bounds_tol) * prob.variable_bounds - x)[:, valid]
    valphas = (valphas / delta[valid]).flatten()

    if prob.bounds.shape[0] > 0:
        # permissible alphas for staying in constraint bounds
        ineqs = prob.inequalities.dot(delta)
        valid = np.abs(ineqs) > sampler.feasibility_tol
        balphas = ((1.0 - sampler.bounds_tol) * prob.bounds - prob.inequalities.dot(x))[
            :, valid
        ]
        balphas = (balphas / ineqs[valid]).flatten()

        # combined alphas
        alphas = np.hstack([valphas, balphas])
    else:
        alphas = valphas

    pos_alphas = alphas[alphas > 0.0]
    neg_alphas = alphas[alphas <= 0.0]
    alpha_range = np.array(
        [
            neg_alphas.max() if len(neg_alphas) > 0 else 0,
            pos_alphas.min() if len(pos_alphas) > 0 else 0,
        ]
    )

    if fraction:
        alpha = alpha_range[0] + fraction * (alpha_range[1] - alpha_range[0])
    else:
        alpha = np.random.uniform(alpha_range[0], alpha_range[1])

    p = x + alpha * delta

    # Numerical instabilities may cause bounds invalidation
    # reset sampler and sample from one of the original warmup directions
    # if that occurs. Also reset if we got stuck.
    if (
        np.any(sampler._bounds_dist(p) < -sampler.bounds_tol)
        or np.abs(np.abs(alpha_range).max() * delta).max() < sampler.bounds_tol
    ):
        if tries > MAX_TRIES:
            raise RuntimeError(
                "Cannot escape sampling region, model seems to be "
                "numerically unstable. Reporting the model to "
                "https://github.com/opencobra/cobrapy/issues "
                "will help us to fix this."
            )
        logger.info("Found bounds infeasibility in sample, resetting to center.")
        newdir = sampler.warmup[np.random.randint(sampler.n_warmup)]
        sampler.retries += 1

        return step(sampler, sampler.center, newdir - sampler.center, None, tries + 1)
    return p
