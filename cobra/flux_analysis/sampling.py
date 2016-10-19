import numpy as np
from ..solvers import solver_dict, get_solver_name


class HRsampler(object):
    """The abstract base class for hit-and-run samplers

    """

    def __init__(self, model, n):
        self.model = model
        self.n_samples = 0
        self.n = n
        self.bounds = np.array([[r.lower_bound, r.upper_bound]
                               for r in model.reactions]).T
        self.samples = np.zeros((self.n, len(self.model.reactions)))

    def generate_warmup(self, n, solver=None, **solver_args):
        if n <= self.n:
            raise ValueError("There must be less warmup points than the \
                              total sample size.")

        solver = solver_dict[get_solver_name() if solver is None else solver]
        lp = solver.create_problem(self.model)
        for i, r in enumerate(self.model.reactions):
            solver.change_variable_objective(lp, i, 0.0)

        rids = np.random.randint(0, len(self.model.reactions), np.ceil(n/2.0))
        for sense in ("minimize", "maximize"):
            for ri in rids:
                solver.change_variable_objective(lp, ri, 1.0)
                solver.solve_problem(lp, objective_sense=sense, **solver_args)
                if self.n_samples == n:
                    break
                self.samples[i, ] = solver.format_solution(lp).x
                self.n_samples += 1
                # revert objective
                solver.change_variable_objective(lp, i, 0.)

    def alpha_bounds(self, x, delta):
        alphas = (self.bounds - x) / delta
        return np.array(alphas[0, ].max(), alphas[1, ].min())

    def sample(n):
        pass
