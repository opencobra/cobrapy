# -*- coding: utf-8 -*-

from __future__ import absolute_import

import optlang.interface


class OptimizationError(Exception):
    def __init__(self, message):
        super(OptimizationError, self).__init__(message)


class Infeasible(OptimizationError):
    pass


class Unbounded(OptimizationError):
    pass


class FeasibleButNotOptimal(OptimizationError):
    pass


class UndefinedSolution(OptimizationError):
    pass


class SolverNotFound(Exception):
    """A simple Exception when a solver can not be found."""

    pass


OPTLANG_TO_EXCEPTIONS_DICT = dict((
    (optlang.interface.INFEASIBLE, Infeasible),
    (optlang.interface.UNBOUNDED, Unbounded),
    (optlang.interface.FEASIBLE, FeasibleButNotOptimal),
    (optlang.interface.UNDEFINED, UndefinedSolution)))
