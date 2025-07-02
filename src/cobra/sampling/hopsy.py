"""Provide hopsy samplers."""

from typing import TYPE_CHECKING, Optional

if TYPE_CHECKING:
    from cobra import Model

from .hr_sampler import HRSampler

import hopsy
import numpy as np
import pandas as pd

class HopsySampler(HRSampler):
    #def __init__(self, model, sampler, processes=processes, thinning=thinning, seed=seed, rounding=True, n_chains=4, **kwargs):
    def __init__(
        self,
        model: "Model",
        sampler,
        thinning: int = 100,
        nproj: Optional[int] = None,
        seed: Optional[int] = None,
        rounding: bool = True,
        processes: int = 1,
        n_chains: int = 4,
        **kwargs,
    ) -> None:
        """Initialize a new HopsySampler."""
        super().__init__(model, thinning, nproj=nproj, seed=seed, **kwargs)

        if self.problem.inequalities.shape[0] > 0:
            A = self.problem.inequalities
            b = self.problem.bounds

            problem = hopsy.Problem(A, b)
            problem = hopsy.add_box_constraints(problem, self.problem.variable_bounds[0], self.problem.variable_bounds[1], simplify=False)
            problem = hopsy.add_equality_constraints(problem, self.problem.equalities, self.problem.b)
            problem = hopsy.round(problem) if rounding else problem
        else:
            # add a dummy inequality constraint which hopefully doesn't intersect with the rest of the polytope
            A = 1e-5*np.ones((1, self.problem.equalities.shape[1]))
            b = 1e+5*np.ones(1)

            problem = hopsy.Problem(A, b)
            problem = hopsy.add_box_constraints(problem, self.problem.variable_bounds[0], self.problem.variable_bounds[1], simplify=False)

            # remove the dummy constraint again
            Ar, br = problem.A[1:], problem.b[1:]
            problem = hopsy.Problem(Ar, br)

            problem = hopsy.add_equality_constraints(problem, self.problem.equalities, self.problem.b)
            problem = hopsy.round(problem) if rounding else problem

        self.processes = processes

        self.mcs = [hopsy.MarkovChain(problem, sampler) for i in range(n_chains)]
        self.rngs = [hopsy.RandomNumberGenerator(self._seed, i) for i in range(n_chains)]

    def sample(self, n: int, fluxes: bool = True) -> pd.DataFrame:
        """
            
        """

        _, samples = hopsy.sample(self.mcs, self.rngs, n_samples=n // len(self.mcs), thinning=self.thinning, n_procs=self.processes, )
        samples = samples.reshape(-1, samples.shape[-1])

        if fluxes:
            names = [r.id for r in self.model.reactions]

            return pd.DataFrame(
                samples[:, self.fwd_idx] - samples[:, self.rev_idx],
                columns=names,
            )
        else:
            names = [v.name for v in self.model.variables]

            return pd.DataFrame(samples, columns=names)

