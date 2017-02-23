# -*- coding: utf-8 -*-

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


class OptimizationError(Exception):
    def __init__(self, message):
        super(OptimizationError, self).__init__(message)


class DefunctError(Exception):
    """Exception for retired functionality

    Parameters
    ----------
    what : string
        The name of the retired object
    alternative : string
        Suggestion for an alternative
    url : string
        A url to alternative resource
    """

    def __init__(self, what, alternative=None, url=None):
        message = "{} has been removed from cobrapy".format(what)
        if alternative is None:
            message += (" without replacement. Raise an issue at "
                        "https://github.com/opencobra/cobrapy if you miss it.")
        if alternative is not None:
            message += ". Consider using '{}' instead".format(alternative)
        if url is not None:
            message += " [{}]".format(url)
        super(DefunctError, self).__init__(message)
