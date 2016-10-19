import numpy as np
from ..solvers import solver_dict, get_solver_name
from copy import deepcopy

BTOL = np.finfo(np.float32).eps
FTOL = BTOL


class HRSampler(object):
    """The abstract base class for hit-and-run samplers

    """

    def __init__(self, model):
        self.model = model
        self.n_samples = 0
        self.bounds = np.array([[r.lower_bound, r.upper_bound]
                               for r in model.reactions]).T

    def generate_fva_warmup(self, solver=None, **solver_args):
        solver = solver_dict[get_solver_name() if solver is None else solver]
        lp = solver.create_problem(self.model)
        for i, r in enumerate(self.model.reactions):
            solver.change_variable_objective(lp, i, 0.0)

        self.n_warmup = 2 * len(self.model.reactions)
        self.warmup = np.zeros((self.n_warmup, len(self.model.reactions)))

        idx = 0
        for sense in ("minimize", "maximize"):
            for i, r in enumerate(self.model.reactions):
                solver.change_variable_objective(lp, i, 1.0)
                solver.solve_problem(lp, objective_sense=sense, **solver_args)
                sol = solver.format_solution(lp, self.model).x
                # some solvers do not enforce bounds too much -> we recontrain
                sol = np.maximum(sol, self.bounds[0, ])
                sol = np.minimum(sol, self.bounds[1, ])
                self.warmup[idx, ] = sol
                idx += 1
                # revert objective
                solver.change_variable_objective(lp, i, 0.)

    def step(self, x, delta):
        delta = delta / np.sqrt((delta * delta).sum())
        nonzero = np.abs(delta) > FTOL
        alphas = ((1.0 - BTOL) * self.bounds - x)[:, nonzero]
        alphas = (alphas / delta[nonzero]).flatten()
        alpha = np.random.uniform(alphas[alphas > 0].min(),
                                  alphas[alphas < 0].max())
        delta[np.logical_not(nonzero)] = 0.0
        return alpha * delta

    def sample(self):
        pass

    def batch(self, batch_size, batch_num):
        for i in range(batch_num):
            b = np.zeros((batch_size, len(self.model.reactions)))
            for j in range(batch_size):
                b[j, ] = self.sample()
            yield b

    def validate(self, samples):
        S = deepcopy(self.model).to_array_based_model().S
        feasibility = np.abs(S.dot(samples.T)).max(axis=0)
        lb_error = (samples - self.bounds[0, ]).min(axis=1)
        ub_error = (self.bounds[1, ] - samples).min(axis=1)

        valid = (feasibility < FTOL) & (lb_error > -BTOL) & (ub_error > -BTOL)
        return valid


class ARCHSampler(HRSampler):
    def __init__(self, model, solver=None, **solver_kwargs):
        super(ARCHSampler, self).__init__(model)
        self.generate_fva_warmup(solver, **solver_kwargs)
        self.prev = self.center = self.warmup.mean(axis=0)

    def sample(self):
        pi = np.random.randint(self.n_warmup)
        # mix in the original warmup points to not get stuck
        p = self.warmup[pi, ]
        delta = p - self.center
        self.prev += self.step(self.prev, delta)
        self.center = (self.n_samples * self.center + self.prev) / (
                       self.n_samples + 1)
        self.n_samples += 1
        return self.prev
