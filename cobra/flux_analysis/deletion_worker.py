# -*- coding: utf-8 -*-

from __future__ import absolute_import

from multiprocessing import Process, Queue, cpu_count

from six import iteritems

from cobra.solvers import get_solver_name, solver_dict


def compute_fba_deletion_worker(cobra_model, solver, job_queue, output_queue,
                                **kwargs):
    solver = solver_dict[get_solver_name() if solver is None else solver]
    lp = solver.create_problem(cobra_model)
    solver_args = kwargs
    solver.solve_problem(lp)
    while True:
        indexes, label = job_queue.get()
        label = indexes if label is None else label
        result = compute_fba_deletion(lp, solver, cobra_model, indexes,
                                      **solver_args)
        output_queue.put((label, result))


def compute_fba_deletion(lp, solver_object, model, indexes, **kwargs):
    s = solver_object
    old_bounds = {}
    for i in indexes:
        reaction = model.reactions[i]
        old_bounds[i] = (reaction.lower_bound, reaction.upper_bound)
        s.change_variable_bounds(lp, i, 0., 0.)
    try:
        s.solve_problem(lp, **kwargs)
    except Exception as e:
        return RuntimeError("solver failure when deleting %s: %s" %
                            (str(indexes), repr(e)))
    status = s.get_status(lp)
    objective = s.get_objective_value(lp) if status == "optimal" else 0.

    # reset the problem, which must be done after reading the solution
    for index, bounds in iteritems(old_bounds):
        s.change_variable_bounds(lp, index, bounds[0], bounds[1])

    if status == "infeasible" or status == "optimal":
        return objective
    else:
        return RuntimeError("solver failure (status %s) for when deleting %s" %
                            (status, str(indexes)))


class CobraDeletionPool(object):
    """A pool of workers for solving deletions

    submit jobs to the pool using submit and recieve results using receive_all
    """
    # Having an existing basis makes solving an existing LP much faster. The
    # most efficient approach is to have a worker function which modifies an LP
    # object and reverts it back after each calculation. Each lp object stores
    # the basis so subsequent LP's are solved more quickely, and memory does
    # not need to be re-allocated each time to create a new problem. Because
    # state is being saved, the workers in the deletion pool are careful about
    # reverting the object after simulating a deletion, and are written to be
    # flexible enough so they can be used in most applications instead of
    # writing a custom worker each time.

    def __init__(self, cobra_model, n_processes=None, solver=None, **kwargs):
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
            p = Process(target=compute_fba_deletion_worker,
                        args=[cobra_model, solver,
                              self.job_queue, self.output_queue],
                        kwargs=kwargs)
            self.processes.append(p)

    def start(self):
        for p in self.processes:
            p.start()

    def terminate(self):
        for p in self.processes:
            p.terminate()

    def __enter__(self):
        self.start()
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        try:
            self.terminate()
        except:
            pass

    def submit(self, indexes, label=None):
        self.job_queue.put((indexes, label))
        self.n_submitted += 1

    def receive_one(self):
        """This function blocks"""
        self.n_complete += 1
        result = self.output_queue.get()
        if isinstance(result[1], Exception):
            raise result[1]
        return result

    def receive_all(self):
        while self.n_complete < self.n_submitted:
            self.n_complete += 1
            result = self.output_queue.get()
            if isinstance(result[1], Exception):
                raise result[1]
            yield result

    @property
    def pids(self):
        return [p.pid for p in self.processes]

    def __del__(self):
        for process in self.processes:
            process.terminate()
            process.join()


class CobraDeletionMockPool(object):
    """Mock pool solves LP's in the same process"""

    def __init__(self, cobra_model, n_processes=1, solver=None, **kwargs):
        if n_processes != 1:
            from warnings import warn
            warn("Mock Pool does not do multiprocessing")
        self.job_queue = []
        self.solver_args = kwargs
        solver_name = get_solver_name() if solver is None else solver
        self.solver = solver_dict[solver_name]
        self.lp = self.solver.create_problem(cobra_model)
        self.solver.solve_problem(self.lp)
        self.model = cobra_model

    def submit(self, indexes, label=None):
        self.job_queue.append((indexes, label))

    def receive_one(self):
        indexes, label = self.job_queue.pop()
        result = compute_fba_deletion(self.lp, self.solver, self.model,
                                      indexes, **self.solver_args)
        if isinstance(result, Exception):
            raise result
        return (label, result)

    def receive_all(self):
        for i in range(len(self.job_queue)):
            indexes, label = self.job_queue.pop()
            result = compute_fba_deletion(self.lp, self.solver, self.model,
                                          indexes, **self.solver_args)
            if isinstance(result, Exception):
                raise result
            yield (label, result)

    def start(self):
        None

    def terminate(self):
        None

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        None
