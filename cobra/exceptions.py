# -*- coding: utf-8 -*-

from __future__ import absolute_import

from collections import defaultdict

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


OPTLANG_TO_EXCEPTIONS_DICT = defaultdict(lambda: OptimizationError)
OPTLANG_TO_EXCEPTIONS_DICT[optlang.interface.INFEASIBLE] = Infeasible
OPTLANG_TO_EXCEPTIONS_DICT[optlang.interface.UNBOUNDED] = Unbounded
OPTLANG_TO_EXCEPTIONS_DICT[optlang.interface.FEASIBLE] = FeasibleButNotOptimal
OPTLANG_TO_EXCEPTIONS_DICT[optlang.interface.UNDEFINED] = UndefinedSolution
