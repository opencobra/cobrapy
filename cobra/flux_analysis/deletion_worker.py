from multiprocessing import Queue, Process, cpu_count

from ..solvers import get_solver_name, solver_dict
from ..external.six import iteritems


def compute_fba_deletion_worker(cobra_model, solver, job_queue, output_queue):
    solver = solver_dict[get_solver_name() if solver is None else solver]
    lp = solver.create_problem(cobra_model)
    solver.solve_problem(lp)
    while True:
        indexes, label = job_queue.get()
        label = indexes if label is None else label
        result = compute_fba_deletion(lp, solver, cobra_model, indexes)
        output_queue.put((label, result))


def compute_fba_deletion(lp, solver_object, model, indexes):
    s = solver_object
    old_bounds = {}
    for i in indexes:
        reaction = model.reactions[i]
        old_bounds[i] = (reaction.lower_bound, reaction.upper_bound)
        s.change_variable_bounds(lp, i, 0., 0.)
    s.solve_problem(lp)
    # reset the problem
    for index, bounds in iteritems(old_bounds):
        s.change_variable_bounds(lp, index, bounds[0], bounds[1])
    return s.get_objective_value(lp) if s.get_status(lp) == "optimal" else 0.


class CobraDeletionPool(object):
    """A pool of workers for solving deletions

    submit jobs to the pool using submit and recieve results using receive_all
    """
    # Having an existing basis makes solving an existing LP much faster. The most
    # efficient approach is to have a worker function which modifies an LP object
    # and reverts it back after each calculation. Each lp object stores the basis
    # so subsequent LP's are solved more quickely, and memory does not need to be
    # re-allocated each time to create a new problem. Because state is being saved,
    # the workers in the deletion pool are careful about reverting the object after
    # simulating a deletion, and are written to be flexible enough so they can be used
    # in most applications instead of writing a custom worker each time.
    def __init__(self, cobra_model, n_processes=None, solver=None):
        if n_processes is None:
            n_processes = min(cpu_count(), 4)
        # start queues
        self.job_queue = Queue()  # format is (indexes, job_label)
        self.n_submitted = 0
        self.n_complete = 0
        self.output_queue = Queue()  # format is (job_label, growth_rate)
        # start processes
        self.processes = []
        for i in range(n_processes):
            p =  Process(target=compute_fba_deletion_worker,
                         args=[cobra_model, solver, self.job_queue, self.output_queue])
            p.start()
            self.processes.append(p)

    def submit(self, indexes, label=None):
        self.job_queue.put((indexes, label))
        self.n_submitted += 1


    def receive_one(self):
        """This function blocks"""
        self.n_complete += 1
        return self.output_queue.get()

    def receive_all(self):
        while self.n_complete < self.n_submitted:
            self.n_complete += 1
            yield self.output_queue.get()

    @property
    def pids(self):
        return [p.pid for p in self.processes]


    def __del__(self):
        for process in self.processes:
            process.terminate()
            process.join()

class CobraDeletionMockPool(object):
    """Mock pool solves LP's in the same process"""
    def __init__(self, cobra_model, n_processes=1, solver=None):
        if n_processes != 1:
            from warnings import warn
            warn("Mock Pool does not do multiprocessing")
        self.job_queue = []
        self.solver = solver_dict[get_solver_name() if solver is None else solver]
        self.lp = self.solver.create_problem(cobra_model)
        self.solver.solve_problem(self.lp)
        self.model = cobra_model

    def submit(self, indexes, label=None):
        self.job_queue.append((indexes, label))

    def receive_one(self):
        indexes, label = self.job_queue.pop()
        return (label, compute_fba_deletion(self.lp, self.solver, self.model, indexes))

    def receive_all(self):
        for i in range(len(self.job_queue)):
            indexes, label = self.job_queue.pop()
            yield (label, compute_fba_deletion(self.lp, self.solver, self.model, indexes))

