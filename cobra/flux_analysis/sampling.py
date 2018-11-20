# -*- coding: utf-8 -*-

"""Module implementing flux sampling for cobra models.

New samplers should derive from the abstract `HRSampler` class
where possible to provide a uniform interface.
"""

from __future__ import absolute_import, division

import ctypes
from collections import namedtuple
from logging import getLogger
from multiprocessing import Array, Pool
from time import time

import numpy as np
import pandas
from optlang.interface import OPTIMAL
from optlang.symbolics import Zero

from cobra.util import (
    constraint_matrices, create_stoichiometric_matrix, nullspace)


LOGGER = getLogger(__name__)
"""The logger for the package."""

bounds_tol = 1e-6
"""The tolerance used for checking bounds feasibility."""

feasibility_tol = 1e-6
"""The tolerance used for checking equalities feasibility."""

max_tries = 100
"""Maximum number of retries for sampling."""

Problem = namedtuple("Problem",
                     ["equalities", "b", "inequalities", "bounds",
                      "variable_fixed", "variable_bounds", "nullspace",
                      "homogeneous"])
"""Defines the matrix representation of a sampling problem.

Attributes
----------
equalities : numpy.array
    All equality constraints in the model.
b : numpy.array
    The right side of the equality constraints.
inequalities : numpy.array
    All inequality constraints in the model.
bounds : numpy.array
    The lower and upper bounds for the inequality constraints.
variable_bounds : numpy.array
    The lower and upper bounds for the variables.
homogeneous: boolean
    Indicates whether the sampling problem is homogenous, e.g. whether there
    exist no non-zero fixed variables or constraints.
nullspace : numpy.matrix
    A matrix containing the nullspace of the equality constraints. Each column
    is one basis vector.
"""


def mp_init(obj):
    """Initialize the multiprocessing pool."""
    global sampler
    sampler = obj


def shared_np_array(shape, data=None, integer=False):
    """Create a new numpy array that resides in shared memory.

    Parameters
    ----------
    shape : tuple of ints
        The shape of the new array.
    data : numpy.array
        Data to copy to the new array. Has to have the same shape.
    integer : boolean
        Whether to use an integer array. Defaults to False which means
        float array.
    """
    size = np.prod(shape)
    if integer:
        array = Array(ctypes.c_int64, int(size))
        np_array = np.frombuffer(array.get_obj(), dtype="int64")
    else:
        array = Array(ctypes.c_double, int(size))
        np_array = np.frombuffer(array.get_obj())
    np_array = np_array.reshape(shape)

    if data is not None:
        if len(shape) != len(data.shape):
            raise ValueError("`data` must have the same dimensions"
                             "as the created array.")
        same = all(x == y for x, y in zip(shape, data.shape))
        if not same:
            raise ValueError("`data` must have the same shape"
                             "as the created array.")
        np_array[:] = data

    return np_array


# Has to be declared outside of class to be used for multiprocessing :(
def _step(sampler, x, delta, fraction=None, tries=0):
    """Sample a new feasible point from the point `x` in direction `delta`."""
    prob = sampler.problem
    valid = ((np.abs(delta) > feasibility_tol) &
             np.logical_not(prob.variable_fixed))
    # permissible alphas for staying in variable bounds
    valphas = ((1.0 - bounds_tol) * prob.variable_bounds - x)[:, valid]
    valphas = (valphas / delta[valid]).flatten()
    if prob.bounds.shape[0] > 0:
        # permissible alphas for staying in constraint bounds
        ineqs = prob.inequalities.dot(delta)
        valid = np.abs(ineqs) > feasibility_tol
        balphas = ((1.0 - bounds_tol) * prob.bounds -
                   prob.inequalities.dot(x))[:, valid]
        balphas = (balphas / ineqs[valid]).flatten()
        # combined alphas
        alphas = np.hstack([valphas, balphas])
    else:
        alphas = valphas
    pos_alphas = alphas[alphas > 0.0]
    neg_alphas = alphas[alphas <= 0.0]
    alpha_range = np.array([neg_alphas.max() if len(neg_alphas) > 0 else 0,
                            pos_alphas.min() if len(pos_alphas) > 0 else 0])

    if fraction:
        alpha = alpha_range[0] + fraction * (alpha_range[1] - alpha_range[0])
    else:
        alpha = np.random.uniform(alpha_range[0], alpha_range[1])
    p = x + alpha * delta

    # Numerical instabilities may cause bounds invalidation
    # reset sampler and sample from one of the original warmup directions
    # if that occurs. Also reset if we got stuck.
    if (np.any(sampler._bounds_dist(p) < -bounds_tol) or
            np.abs(np.abs(alpha_range).max() * delta).max() < bounds_tol):
        if tries > max_tries:
            raise RuntimeError("Can not escape sampling region, model seems"
                               " numerically unstable :( Reporting the "
                               "model to "
                               "https://github.com/opencobra/cobrapy/issues "
                               "will help us to fix this :)")
        LOGGER.info("found bounds infeasibility in sample, "
                    "resetting to center")
        newdir = sampler.warmup[np.random.randint(sampler.n_warmup)]
        sampler.retries += 1
        return _step(sampler, sampler.center, newdir - sampler.center, None,
                     tries + 1)
    return p


class HRSampler(object):
    """The abstract base class for hit-and-run samplers.

    Parameters
    ----------
    model : cobra.Model
        The cobra model from which to generate samples.
    thinning : int
        The thinning factor of the generated sampling chain. A thinning of 10
        means samples are returned every 10 steps.
    nproj : int > 0, optional
        How often to reproject the sampling point into the feasibility space.
        Avoids numerical issues at the cost of lower sampling. If you observe
        many equality constraint violations with `sampler.validate` you should
        lower this number.
    seed : int > 0, optional
        The random number seed that should be used.

    Attributes
    ----------
    model : cobra.Model
        The cobra model from which the samples get generated.
    thinning : int
        The currently used thinning factor.
    n_samples : int
        The total number of samples that have been generated by this
        sampler instance.
    retries : int
        The overall of sampling retries the sampler has observed. Larger
        values indicate numerical instabilities.
    problem : collections.namedtuple
        A python object whose attributes define the entire sampling problem in
        matrix form. See docstring of `Problem`.
    warmup : a numpy matrix
        A matrix of with as many columns as reactions in the model and more
        than 3 rows containing a warmup sample in each row. None if no warmup
        points have been generated yet.
    nproj : int
        How often to reproject the sampling point into the feasibility space.
    seed : positive integer, optional
        Sets the random number seed. Initialized to the current time stamp if
        None.
    fwd_idx : np.array
        Has one entry for each reaction in the model containing the index of
        the respective forward variable.
    rev_idx : np.array
        Has one entry for each reaction in the model containing the index of
        the respective reverse variable.
    """

    def __init__(self, model, thinning, nproj=None, seed=None):
        """Initialize a new sampler object."""
        # This currently has to be done to reset the solver basis which is
        # required to get deterministic warmup point generation
        # (in turn required for a working `seed` arg)
        if model.solver.is_integer:
            raise TypeError("sampling does not work with integer problems :(")
        self.model = model.copy()
        self.thinning = thinning
        if nproj is None:
            self.nproj = int(min(len(self.model.variables)**3, 1e6))
        else:
            self.nproj = nproj
        self.n_samples = 0
        self.retries = 0
        self.problem = self.__build_problem()
        # Set up a map from reaction -> forward/reverse variable
        var_idx = {v: idx for idx, v in enumerate(model.variables)}
        self.fwd_idx = np.array([var_idx[r.forward_variable]
                                 for r in model.reactions])
        self.rev_idx = np.array([var_idx[r.reverse_variable]
                                 for r in model.reactions])
        self.warmup = None
        if seed is None:
            self._seed = int(time())
        else:
            self._seed = seed
        # Avoid overflow
        self._seed = self._seed % np.iinfo(np.int32).max

    def __build_problem(self):
        """Build the matrix representation of the sampling problem."""
        # Set up the mathematical problem
        prob = constraint_matrices(self.model, zero_tol=feasibility_tol)
        # check if there any non-zero equality constraints
        equalities = prob.equalities
        b = prob.b
        bounds = np.atleast_2d(prob.bounds).T
        var_bounds = np.atleast_2d(prob.variable_bounds).T
        homogeneous = all(np.abs(b) < feasibility_tol)
        fixed_non_zero = np.abs(prob.variable_bounds[:, 1]) > feasibility_tol
        fixed_non_zero &= prob.variable_fixed
        # check if there are any non-zero fixed variables, add them as
        # equalities to the stoichiometric matrix
        if any(fixed_non_zero):
            n_fixed = fixed_non_zero.sum()
            rows = np.zeros((n_fixed, prob.equalities.shape[1]))
            rows[range(n_fixed), np.where(fixed_non_zero)] = 1.0
            equalities = np.vstack([equalities, rows])
            var_b = prob.variable_bounds[:, 1]
            b = np.hstack([b, var_b[fixed_non_zero]])
            homogeneous = False
        # Set up a projection that can cast point into the nullspace
        nulls = nullspace(equalities)
        # convert bounds to a matrix and add variable bounds as well
        return Problem(
            equalities=shared_np_array(equalities.shape, equalities),
            b=shared_np_array(b.shape, b),
            inequalities=shared_np_array(prob.inequalities.shape,
                                         prob.inequalities),
            bounds=shared_np_array(bounds.shape, bounds),
            variable_fixed=shared_np_array(prob.variable_fixed.shape,
                                           prob.variable_fixed, integer=True),
            variable_bounds=shared_np_array(var_bounds.shape, var_bounds),
            nullspace=shared_np_array(nulls.shape, nulls),
            homogeneous=homogeneous
        )

    def generate_fva_warmup(self):
        """Generate the warmup points for the sampler.

        Generates warmup points by setting each flux as the sole objective
        and minimizing/maximizing it. Also caches the projection of the
        warmup points into the nullspace for non-homogeneous problems (only
        if necessary).
        """
        self.n_warmup = 0
        reactions = self.model.reactions
        self.warmup = np.zeros((2 * len(reactions), len(self.model.variables)))
        self.model.objective = Zero
        for sense in ("min", "max"):
            self.model.objective_direction = sense
            for i, r in enumerate(reactions):
                variables = (self.model.variables[self.fwd_idx[i]],
                             self.model.variables[self.rev_idx[i]])
                # Omit fixed reactions if they are non-homogeneous
                if r.upper_bound - r.lower_bound < bounds_tol:
                    LOGGER.info("skipping fixed reaction %s" % r.id)
                    continue
                self.model.objective.set_linear_coefficients(
                    {variables[0]: 1, variables[1]: -1})
                self.model.slim_optimize()
                if not self.model.solver.status == OPTIMAL:
                    LOGGER.info("can not maximize reaction %s, skipping it" %
                                r.id)
                    continue
                primals = self.model.solver.primal_values
                sol = [primals[v.name] for v in self.model.variables]
                self.warmup[self.n_warmup, ] = sol
                self.n_warmup += 1
                # Reset objective
                self.model.objective.set_linear_coefficients(
                    {variables[0]: 0, variables[1]: 0})
        # Shrink to measure
        self.warmup = self.warmup[0:self.n_warmup, :]
        # Remove redundant search directions
        keep = np.logical_not(self._is_redundant(self.warmup))
        self.warmup = self.warmup[keep, :]
        self.n_warmup = self.warmup.shape[0]

        # Catch some special cases
        if len(self.warmup.shape) == 1 or self.warmup.shape[0] == 1:
            raise ValueError("Your flux cone consists only of a single point!")
        elif self.n_warmup == 2:
            if not self.problem.homogeneous:
                raise ValueError("Can not sample from an inhomogenous problem"
                                 " with only 2 search directions :(")
            LOGGER.info("All search directions on a line, adding another one.")
            newdir = self.warmup.T.dot([0.25, 0.25])
            self.warmup = np.vstack([self.warmup, newdir])
            self.n_warmup += 1

        # Shrink warmup points to measure
        self.warmup = shared_np_array(
            (self.n_warmup, len(self.model.variables)), self.warmup)

    def _reproject(self, p):
        """Reproject a point into the feasibility region.

        This function is guaranteed to return a new feasible point. However,
        no guarantees in terms of proximity to the original point can be made.

        Parameters
        ----------
        p : numpy.array
            The current sample point.

        Returns
        -------
        numpy.array
            A new feasible point. If `p` was feasible it wil return p.
        """
        nulls = self.problem.nullspace
        equalities = self.problem.equalities
        # don't reproject if point is feasible
        if np.allclose(equalities.dot(p), self.problem.b,
                       rtol=0, atol=feasibility_tol):
            new = p
        else:
            LOGGER.info("feasibility violated in sample"
                        " %d, trying to reproject" % self.n_samples)
            new = nulls.dot(nulls.T.dot(p))
        # Projections may violate bounds
        # set to random point in space in that case
        if any(new != p):
            LOGGER.info("reprojection failed in sample"
                        " %d, using random point in space" % self.n_samples)
            new = self._random_point()
        return new

    def _random_point(self):
        """Find an approximately random point in the flux cone."""
        idx = np.random.randint(self.n_warmup,
                                size=min(2, np.ceil(np.sqrt(self.n_warmup))))
        return self.warmup[idx, :].mean(axis=0)

    def _is_redundant(self, matrix, cutoff=1.0 - feasibility_tol):
        """Identify rdeundant rows in a matrix that can be removed."""
        # Avoid zero variances
        extra_col = matrix[:, 0] + 1
        # Avoid zero rows being correlated with constant rows
        extra_col[matrix.sum(axis=1) == 0] = 2
        corr = np.corrcoef(np.c_[matrix, extra_col])
        corr = np.tril(corr, -1)
        return (np.abs(corr) > cutoff).any(axis=1)

    def _bounds_dist(self, p):
        """Get the lower and upper bound distances. Negative is bad."""
        prob = self.problem
        lb_dist = (p - prob.variable_bounds[0, ]).min()
        ub_dist = (prob.variable_bounds[1, ] - p).min()
        if prob.bounds.shape[0] > 0:
            const = prob.inequalities.dot(p)
            const_lb_dist = (const - prob.bounds[0, ]).min()
            const_ub_dist = (prob.bounds[1, ] - const).min()
            lb_dist = min(lb_dist, const_lb_dist)
            ub_dist = min(ub_dist, const_ub_dist)
        return np.array([lb_dist, ub_dist])

    def sample(self, n, fluxes=True):
        """Abstract sampling function.

        Should be overwritten by child classes.
        """
        pass

    def batch(self, batch_size, batch_num, fluxes=True):
        """Create a batch generator.

        This is useful to generate n batches of m samples each.

        Parameters
        ----------
        batch_size : int
            The number of samples contained in each batch (m).
        batch_num : int
            The number of batches in the generator (n).
        fluxes : boolean
            Whether to return fluxes or the internal solver variables. If set
            to False will return a variable for each forward and backward flux
            as well as all additional variables you might have defined in the
            model.

        Yields
        ------
        pandas.DataFrame
            A DataFrame with dimensions (batch_size x n_r) containing
            a valid flux sample for a total of n_r reactions (or variables if
            fluxes=False) in each row.
        """
        for i in range(batch_num):
            yield self.sample(batch_size, fluxes=fluxes)

    def validate(self, samples):
        """Validate a set of samples for equality and inequality feasibility.

        Can be used to check whether the generated samples and warmup points
        are feasible.

        Parameters
        ----------
        samples : numpy.matrix
            Must be of dimension (n_samples x n_reactions). Contains the
            samples to be validated. Samples must be from fluxes.

        Returns
        -------
        numpy.array
            A one-dimensional numpy array of length containing
            a code of 1 to 3 letters denoting the validation result:

            - 'v' means feasible in bounds and equality constraints
            - 'l' means a lower bound violation
            - 'u' means a lower bound validation
            - 'e' means and equality constraint violation
        """
        samples = np.atleast_2d(samples)
        prob = self.problem

        if samples.shape[1] == len(self.model.reactions):
            S = create_stoichiometric_matrix(self.model)
            b = np.array([self.model.constraints[m.id].lb for m in
                          self.model.metabolites])
            bounds = np.array([r.bounds for r in self.model.reactions]).T
        elif samples.shape[1] == len(self.model.variables):
            S = prob.equalities
            b = prob.b
            bounds = prob.variable_bounds
        else:
            raise ValueError("Wrong number of columns. samples must have a "
                             "column for each flux or variable defined in the "
                             "model!")

        feasibility = np.abs(S.dot(samples.T).T - b)
        feasibility = feasibility.max(axis=1)
        lb_error = (samples - bounds[0, ]).min(axis=1)
        ub_error = (bounds[1, ] - samples).min(axis=1)

        if (samples.shape[1] == len(self.model.variables) and
                prob.inequalities.shape[0]):
            consts = prob.inequalities.dot(samples.T)
            lb_error = np.minimum(
                lb_error,
                (consts - prob.bounds[0, ]).min(axis=1))
            ub_error = np.minimum(
                ub_error,
                (prob.bounds[1, ] - consts).min(axis=1)
            )

        valid = (
            (feasibility < feasibility_tol) &
            (lb_error > -bounds_tol) &
            (ub_error > -bounds_tol))
        codes = np.repeat("", valid.shape[0]).astype(np.dtype((str, 3)))
        codes[valid] = "v"
        codes[lb_error <= -bounds_tol] = np.char.add(
            codes[lb_error <= -bounds_tol], "l")
        codes[ub_error <= -bounds_tol] = np.char.add(
            codes[ub_error <= -bounds_tol], "u")
        codes[feasibility > feasibility_tol] = np.char.add(
            codes[feasibility > feasibility_tol], "e")
        return codes


class ACHRSampler(HRSampler):
    """Artificial Centering Hit-and-Run sampler.

    A sampler with low memory footprint and good convergence.

    Parameters
    ----------
    model : a cobra model
        The cobra model from which to generate samples.
    thinning : int, optional
        The thinning factor of the generated sampling chain. A thinning of 10
        means samples are returned every 10 steps.
    nproj : int > 0, optional
        How often to reproject the sampling point into the feasibility space.
        Avoids numerical issues at the cost of lower sampling. If you observe
        many equality constraint violations with `sampler.validate` you should
        lower this number.
    seed : int > 0, optional
        Sets the random number seed. Initialized to the current time stamp if
        None.

    Attributes
    ----------
    model : cobra.Model
        The cobra model from which the samples get generated.
    thinning : int
        The currently used thinning factor.
    n_samples : int
        The total number of samples that have been generated by this
        sampler instance.
    problem : collections.namedtuple
        A python object whose attributes define the entire sampling problem in
        matrix form. See docstring of `Problem`.
    warmup : a numpy matrix
        A matrix of with as many columns as reactions in the model and more
        than 3 rows containing a warmup sample in each row. None if no warmup
        points have been generated yet.
    retries : int
        The overall of sampling retries the sampler has observed. Larger
        values indicate numerical instabilities.
    seed : positive integer, optional
        Sets the random number seed. Initialized to the current time stamp if
        None.
    nproj : int
        How often to reproject the sampling point into the feasibility space.
    fwd_idx : np.array
        Has one entry for each reaction in the model containing the index of
        the respective forward variable.
    rev_idx : np.array
        Has one entry for each reaction in the model containing the index of
        the respective reverse variable.
    prev : numpy array
        The current/last flux sample generated.
    center : numpy array
        The center of the sampling space as estimated by the mean of all
        previously generated samples.

    Notes
    -----
    ACHR generates samples by choosing new directions from the sampling space's
    center and the warmup points. The implementation used here is the same
    as in the Matlab Cobra Toolbox [2]_ and uses only the initial warmup points
    to generate new directions and not any other previous iterates. This
    usually gives better mixing since the startup points are chosen to span
    the space in a wide manner. This also makes the generated sampling chain
    quasi-markovian since the center converges rapidly.

    Memory usage is roughly in the order of (2 * number reactions)^2
    due to the required nullspace matrices and warmup points. So large
    models easily take up a few GB of RAM.

    References
    ----------
    .. [1] Direction Choice for Accelerated Convergence in Hit-and-Run Sampling
       David E. Kaufman Robert L. Smith
       Operations Research 199846:1 , 84-95
       https://doi.org/10.1287/opre.46.1.84
    .. [2] https://github.com/opencobra/cobratoolbox
    """

    def __init__(self, model, thinning=100, nproj=None, seed=None):
        """Initialize a new ACHRSampler."""
        super(ACHRSampler, self).__init__(model, thinning, nproj=nproj,
                                          seed=seed)
        self.generate_fva_warmup()
        self.prev = self.center = self.warmup.mean(axis=0)
        np.random.seed(self._seed)

    def __single_iteration(self):
        pi = np.random.randint(self.n_warmup)
        # mix in the original warmup points to not get stuck
        delta = self.warmup[pi, ] - self.center
        self.prev = _step(self, self.prev, delta)
        if self.problem.homogeneous and (self.n_samples *
                                         self.thinning % self.nproj == 0):
            self.prev = self._reproject(self.prev)
            self.center = self._reproject(self.center)
        self.center = ((self.n_samples * self.center) / (self.n_samples + 1) +
                       self.prev / (self.n_samples + 1))
        self.n_samples += 1

    def sample(self, n, fluxes=True):
        """Generate a set of samples.

        This is the basic sampling function for all hit-and-run samplers.

        Parameters
        ----------
        n : int
            The number of samples that are generated at once.
        fluxes : boolean
            Whether to return fluxes or the internal solver variables. If set
            to False will return a variable for each forward and backward flux
            as well as all additional variables you might have defined in the
            model.

        Returns
        -------
        numpy.matrix
            Returns a matrix with `n` rows, each containing a flux sample.

        Notes
        -----
        Performance of this function linearly depends on the number
        of reactions in your model and the thinning factor.
        """
        samples = np.zeros((n, self.warmup.shape[1]))
        for i in range(1, self.thinning * n + 1):
            self.__single_iteration()
            if i % self.thinning == 0:
                samples[i//self.thinning - 1, ] = self.prev

        if fluxes:
            names = [r.id for r in self.model.reactions]
            return pandas.DataFrame(
                samples[:, self.fwd_idx] - samples[:, self.rev_idx],
                columns=names)
        else:
            names = [v.name for v in self.model.variables]
            return pandas.DataFrame(samples, columns=names)


# Unfortunately this has to be outside the class to be usable with
# multiprocessing :()
def _sample_chain(args):
    """Sample a single chain for OptGPSampler.

    center and n_samples are updated locally and forgotten afterwards.
    """
    n, idx = args       # has to be this way to work in Python 2.7
    center = sampler.center
    np.random.seed((sampler._seed + idx) % np.iinfo(np.int32).max)
    pi = np.random.randint(sampler.n_warmup)
    prev = sampler.warmup[pi, ]
    prev = _step(sampler, center, prev - center, 0.95)
    n_samples = max(sampler.n_samples, 1)
    samples = np.zeros((n, center.shape[0]))

    for i in range(1, sampler.thinning * n + 1):
        pi = np.random.randint(sampler.n_warmup)
        delta = sampler.warmup[pi, ] - center

        prev = _step(sampler, prev, delta)
        if sampler.problem.homogeneous and (
                n_samples * sampler.thinning % sampler.nproj == 0):
            prev = sampler._reproject(prev)
            center = sampler._reproject(center)
        if i % sampler.thinning == 0:
            samples[i//sampler.thinning - 1, ] = prev
        center = ((n_samples * center) / (n_samples + 1) +
                  prev / (n_samples + 1))
        n_samples += 1

    return (sampler.retries, samples)


class OptGPSampler(HRSampler):
    """A parallel optimized sampler.

    A parallel sampler with fast convergence and parallel execution. See [1]_
    for details.

    Parameters
    ----------
    model : cobra.Model
        The cobra model from which to generate samples.
    processes: int
        The number of processes used during sampling.
    thinning : int, optional
        The thinning factor of the generated sampling chain. A thinning of 10
        means samples are returned every 10 steps.
    nproj : int > 0, optional
        How often to reproject the sampling point into the feasibility space.
        Avoids numerical issues at the cost of lower sampling. If you observe
        many equality constraint violations with `sampler.validate` you should
        lower this number.
    seed : int > 0, optional
        Sets the random number seed. Initialized to the current time stamp if
        None.

    Attributes
    ----------
    model : cobra.Model
        The cobra model from which the samples get generated.
    thinning : int
        The currently used thinning factor.
    n_samples : int
        The total number of samples that have been generated by this
        sampler instance.
    problem : collections.namedtuple
        A python object whose attributes define the entire sampling problem in
        matrix form. See docstring of `Problem`.
    warmup : a numpy matrix
        A matrix of with as many columns as reactions in the model and more
        than 3 rows containing a warmup sample in each row. None if no warmup
        points have been generated yet.
    retries : int
        The overall of sampling retries the sampler has observed. Larger
        values indicate numerical instabilities.
    seed : positive integer, optional
        Sets the random number seed. Initialized to the current time stamp if
        None.
    nproj : int
        How often to reproject the sampling point into the feasibility space.
    fwd_idx : np.array
        Has one entry for each reaction in the model containing the index of
        the respective forward variable.
    rev_idx : np.array
        Has one entry for each reaction in the model containing the index of
        the respective reverse variable.
    prev : numpy.array
        The current/last flux sample generated.
    center : numpy.array
        The center of the sampling space as estimated by the mean of all
        previously generated samples.

    Notes
    -----
    The sampler is very similar to artificial centering where each process
    samples its own chain. Initial points are chosen randomly from the warmup
    points followed by a linear transformation that pulls the points towards
    the a little bit towards the center of the sampling space.

    If the number of processes used is larger than one the requested
    number of samples is adjusted to the smallest multiple of the number of
    processes larger than the requested sample number. For instance, if you
    have 3 processes and request 8 samples you will receive 9.

    Memory usage is roughly in the order of (2 * number reactions)^2
    due to the required nullspace matrices and warmup points. So large
    models easily take up a few GB of RAM. However, most of the large matrices
    are kept in shared memory. So the RAM usage is independent of the number
    of processes.

    References
    ----------
    .. [1] Megchelenbrink W, Huynen M, Marchiori E (2014)
       optGpSampler: An Improved Tool for Uniformly Sampling the Solution-Space
       of Genome-Scale Metabolic Networks.
       PLoS ONE 9(2): e86587.
       https://doi.org/10.1371/journal.pone.0086587
    """

    def __init__(self, model, processes, thinning=100, nproj=None,
                 seed=None):
        """Initialize a new OptGPSampler."""
        super(OptGPSampler, self).__init__(model, thinning, seed=seed)
        self.generate_fva_warmup()
        self.processes = processes

        # This maps our saved center into shared memory,
        # meaning they are synchronized across processes
        self.center = shared_np_array((len(model.variables), ),
                                      self.warmup.mean(axis=0))

    def sample(self, n, fluxes=True):
        """Generate a set of samples.

        This is the basic sampling function for all hit-and-run samplers.

        Paramters
        ---------
        n : int
            The minimum number of samples that are generated at once
            (see Notes).
        fluxes : boolean
            Whether to return fluxes or the internal solver variables. If set
            to False will return a variable for each forward and backward flux
            as well as all additional variables you might have defined in the
            model.

        Returns
        -------
        numpy.matrix
            Returns a matrix with `n` rows, each containing a flux sample.

        Notes
        -----
        Performance of this function linearly depends on the number
        of reactions in your model and the thinning factor.

        If the number of processes is larger than one, computation is split
        across as the CPUs of your machine. This may shorten computation time.
        However, there is also overhead in setting up parallel computation so
        we recommend to calculate large numbers of samples at once
        (`n` > 1000).
        """
        if self.processes > 1:
            n_process = np.ceil(n / self.processes).astype(int)
            n = n_process * self.processes
            # The cast to list is weird but not doing it gives recursion
            # limit errors, something weird going on with multiprocessing
            args = list(zip(
                [n_process] * self.processes, range(self.processes)))
            # No with statement or starmap here since Python 2.x
            # does not support it :(
            mp = Pool(self.processes, initializer=mp_init, initargs=(self,))
            results = mp.map(_sample_chain, args, chunksize=1)
            mp.close()
            mp.join()
            chains = np.vstack([r[1] for r in results])
            self.retries += sum(r[0] for r in results)
        else:
            mp_init(self)
            results = _sample_chain((n, 0))
            chains = results[1]

        # Update the global center
        self.center = (self.n_samples * self.center +
                       np.atleast_2d(chains).sum(0)) / (self.n_samples + n)
        self.n_samples += n

        if fluxes:
            names = [r.id for r in self.model.reactions]
            return pandas.DataFrame(
                chains[:, self.fwd_idx] - chains[:, self.rev_idx],
                columns=names)
        else:
            names = [v.name for v in self.model.variables]
            return pandas.DataFrame(chains, columns=names)

    # Models can be large so don't pass them around during multiprocessing
    def __getstate__(self):
        """Return the object for serialization."""
        d = dict(self.__dict__)
        del d['model']
        return d


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
    seed : positive integer, optional
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
