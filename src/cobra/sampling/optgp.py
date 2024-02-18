"""Provide the OptGP sampler class and helper functions."""

from typing import TYPE_CHECKING, Dict, Optional, Tuple

import numpy as np
import pandas as pd

from ..core.configuration import Configuration
from ..util import ProcessPool
from .core import step
from .hr_sampler import HRSampler, shared_np_array


if TYPE_CHECKING:
    from cobra import Model


__all__ = ("OptGPSampler",)


configuration = Configuration()


class OptGPSampler(HRSampler):
    """
    Improved Artificial Centering Hit-and-Run sampler.

    A parallel sampler with fast convergence and parallel execution.
    See [1]_ for details.

    Parameters
    ----------
    model : cobra.Model
        The cobra model from which to generate samples.
    processes: int, optional
        The number of processes used during sampling
        (default cobra.Configuration.processes).
    thinning : int, optional
        The thinning factor of the generated sampling chain. A thinning of
        10 means samples are returned every 10 steps (default 100).
    nproj : int > 0, optional
        How often to reproject the sampling point into the feasibility
        space. Avoids numerical issues at the cost of lower sampling. If
        you observe many equality constraint violations with
        `sampler.validate` you should lower this number (default None).
    seed : int > 0, optional
        Sets the random number seed. Initialized to the current time stamp
        if None (default None).

    Attributes
    ----------
    n_samples : int
        The total number of samples that have been generated by this
        sampler instance.
    problem : typing.NamedTuple
        A NamedTuple whose attributes define the entire sampling problem in
        matrix form.
    warmup : numpy.matrix
        A numpy matrix with as many columns as reactions in the model and
        more than 3 rows containing a warmup sample in each row. None if no
        warmup points have been generated yet.
    retries : int
        The overall of sampling retries the sampler has observed. Larger
        values indicate numerical instabilities.
    fwd_idx : numpy.array
        A numpy array having one entry for each reaction in the model,
        containing the index of the respective forward variable.
    rev_idx : numpy.array
        A numpy array having one entry for each reaction in the model,
        containing the index of the respective reverse variable.
    prev : numpy.array
        The current/last flux sample generated.
    center : numpy.array
        The center of the sampling space as estimated by the mean of all
        previously generated samples.

    Notes
    -----
    The sampler is very similar to artificial centering where each process
    samples its own chain. Initial points are chosen randomly from the
    warmup points followed by a linear transformation that pulls the points
    a little bit towards the center of the sampling space.

    If the number of processes used is larger than the one requested,
    number of samples is adjusted to the smallest multiple of the number of
    processes larger than the requested sample number. For instance, if you
    have 3 processes and request 8 samples, you will receive 9.

    Memory usage is roughly in the order of (2 * number of reactions)^2
    due to the required nullspace matrices and warmup points. So, large
    models easily take up a few GBs of RAM. However, most of the large
    matrices are kept in shared memory. So the RAM usage is independent of
    the number of processes.

    References
    ----------
    .. [1] Megchelenbrink W, Huynen M, Marchiori E (2014)
       optGpSampler: An Improved Tool for Uniformly Sampling the Solution-Space
       of Genome-Scale Metabolic Networks.
       PLoS ONE 9(2): e86587.
       https://doi.org/10.1371/journal.pone.0086587

    """

    def __init__(
        self,
        model: "Model",
        thinning: int = 100,
        processes: Optional[int] = None,
        nproj: Optional[int] = None,
        seed: Optional[int] = None,
        **kwargs,
    ) -> None:
        """Initialize a new OptGPSampler."""
        super().__init__(model, thinning, nproj=nproj, seed=seed, **kwargs)
        self.generate_fva_warmup()

        if processes is None:
            self.processes = configuration.processes
        else:
            self.processes = processes

        # This maps our saved center into shared memory,
        # meaning they are synchronized across processes
        self.center = shared_np_array(
            (len(self.model.variables),), self.warmup.mean(axis=0)
        )

    def sample(self, n: int, fluxes: bool = True) -> pd.DataFrame:
        """Generate a set of samples.

        This is the basic sampling function for all hit-and-run samplers.

        Parameters
        ----------
        n : int
            The minimum number of samples that are generated at once.
        fluxes : bool, optional
            Whether to return fluxes or the internal solver variables. If
            set to False, will return a variable for each forward and
            backward flux as well as all additional variables you might
            have defined in the model (default True).

        Returns
        -------
        pandas.DataFrame
            Returns a pandas DataFrame with `n` rows, each containing a
            flux sample.

        Notes
        -----
        Performance of this function linearly depends on the number
        of reactions in your model and the thinning factor.

        If the number of processes is larger than one, computation is split
        across the CPU cores of your machine. This may shorten computation
        time. However, there is also overhead in setting up parallel
        computation primitives so, we recommend to calculate large numbers
        of samples at once (`n` > 1000).

        """
        if self.processes > 1:
            n_process = np.ceil(n / self.processes).astype(int)
            n = n_process * self.processes

            # The cast to list is weird but not doing it gives recursion
            # limit errors, something weird going on with multiprocessing
            args = list(zip([n_process] * self.processes, range(self.processes)))

            with ProcessPool(
                self.processes, initializer=mp_init, initargs=(self,)
            ) as pool:
                results = pool.map(_sample_chain, args, chunksize=1)

            chains = np.vstack([r[1] for r in results])
            self.retries += sum(r[0] for r in results)
        else:
            mp_init(self)
            results = _sample_chain((n, 0))
            chains = results[1]

        # Update the global center
        self.center = (self.n_samples * self.center + np.atleast_2d(chains).sum(0)) / (
            self.n_samples + n
        )
        self.n_samples += n

        if fluxes:
            names = [r.id for r in self.model.reactions]

            return pd.DataFrame(
                chains[:, self.fwd_idx] - chains[:, self.rev_idx],
                columns=names,
            )
        else:
            names = [v.name for v in self.model.variables]

            return pd.DataFrame(chains, columns=names)

    # Models can be large so don't pass them around during multiprocessing
    def __getstate__(self) -> Dict:
        """Return the object for serialization."""
        d = dict(self.__dict__)
        del d["model"]
        return d


def mp_init(obj: OptGPSampler) -> None:
    """Initialize the multiprocessing pool."""
    global sampler
    sampler = obj


# Unfortunately this has to be outside the class to be usable with
# multiprocessing :()
def _sample_chain(args: Tuple[int, int]) -> Tuple[int, OptGPSampler]:
    """Sample a single chain for OptGPSampler.

    `center` and `n_samples` are updated locally and forgotten afterwards.

    """
    n, idx = args  # has to be this way to work in Python 2.7
    center = sampler.center
    np.random.seed((sampler._seed + idx) % np.iinfo(np.int32).max)
    pi = np.random.randint(sampler.n_warmup)

    prev = sampler.warmup[pi, :]
    prev = step(sampler, center, prev - center, 0.95)

    n_samples = max(sampler.n_samples, 1)
    samples = np.zeros((n, center.shape[0]))

    for i in range(1, sampler.thinning * n + 1):
        pi = np.random.randint(sampler.n_warmup)
        delta = sampler.warmup[pi, :] - center

        prev = step(sampler, prev, delta)

        if sampler.problem.homogeneous and (
            n_samples * sampler.thinning % sampler.nproj == 0
        ):
            prev = sampler._reproject(prev)
            center = sampler._reproject(center)

        if i % sampler.thinning == 0:
            samples[i // sampler.thinning - 1, :] = prev

        center = (n_samples * center) / (n_samples + 1) + prev / (n_samples + 1)
        n_samples += 1

    return (sampler.retries, samples)
