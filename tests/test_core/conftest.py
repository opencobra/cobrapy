"""Define module level fixtures."""

from typing import TYPE_CHECKING, Tuple

import pytest

from cobra.util.solver import solvers


if TYPE_CHECKING:
    from cobra import Model, Solution


solver_trials = [
    "glpk",
    pytest.param(
        "cplex",
        marks=pytest.mark.skipif(
            "cplex" not in solvers, reason="No CPLEX found on PYTHONPATH"
        ),
    ),
    pytest.param(
        "gurobi",
        marks=pytest.mark.skipif(
            "gurobi" not in solvers, reason="No Gurobi found on PYTHONPATH"
        ),
    ),
]


@pytest.fixture(scope="function", params=solver_trials)
def solved_model(
    request: pytest.FixtureRequest, model: "Model"
) -> Tuple["Solution", "Model"]:
    """Return solved model."""
    model.solver = request.param
    solution = model.optimize()
    return solution, model
