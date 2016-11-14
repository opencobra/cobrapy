from __future__ import absolute_import, print_function

import optlang.interface


class SolveError(Exception):
    def __init__(self, message):
        super(SolveError, self).__init__(message)


class Infeasible(SolveError):
    pass


class Unbounded(SolveError):
    pass


class FeasibleButNotOptimal(SolveError):
    pass


class UndefinedSolution(SolveError):
    pass


_OPTLANG_TO_EXCEPTIONS_DICT = dict((
    (optlang.interface.INFEASIBLE, Infeasible),
    (optlang.interface.UNBOUNDED, Unbounded),
    (optlang.interface.FEASIBLE, FeasibleButNotOptimal),
    (optlang.interface.UNDEFINED, UndefinedSolution)))
