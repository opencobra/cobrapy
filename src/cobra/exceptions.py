"""Module for shared exceptions in the Cobra package."""
import optlang.interface


class OptimizationError(Exception):
    """Exception for Optimization issues."""

    def __init__(self, message):
        """Inherit parent behaviors."""
        super(OptimizationError, self).__init__(message)


class Infeasible(OptimizationError):
    """Exception for Infeasible issues."""

    pass


class Unbounded(OptimizationError):
    """Exception for Unbounded issues."""

    pass


class FeasibleButNotOptimal(OptimizationError):
    """Exception for Non-Optimal issues."""

    pass


class UndefinedSolution(OptimizationError):
    """Exception for Undefined issues."""

    pass


class SolverNotFound(Exception):
    """A simple Exception when a solver can not be found."""

    pass


OPTLANG_TO_EXCEPTIONS_DICT = dict(
    (
        (optlang.interface.INFEASIBLE, Infeasible),
        (optlang.interface.UNBOUNDED, Unbounded),
        (optlang.interface.FEASIBLE, FeasibleButNotOptimal),
        (optlang.interface.UNDEFINED, UndefinedSolution),
    )
)
