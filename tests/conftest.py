"""Define global fixtures."""

from pathlib import Path
from pickle import load as _load
from typing import List, Tuple

import importlib_resources
import pytest

import cobra
from cobra import Metabolite, Model, Reaction, Solution
from cobra.io import read_sbml_model
from cobra.util import solver as sutil


data_dir = Path(__file__).parent / "data"
cobra_data_dir = importlib_resources.files(cobra.data)


def create_test_model(model_name: str = "salmonella") -> Model:
    """Return a cobra model for testing.

    Parameters
    ----------
    model_name: str
        One of 'ecoli', 'textbook', or 'salmonella', or the
        path to a pickled cobra.Model .

    Returns
    -------
    cobra.Model
        The cobra model.

    Raises
    ------
    OSError
        If no file is found at the path specified by `model_name`.

    """
    if model_name == "ecoli":
        ecoli_sbml = str((cobra_data_dir / "iJO1366.xml.gz").resolve())
        return read_sbml_model(ecoli_sbml)
    elif model_name == "textbook":
        textbook_sbml = str((cobra_data_dir / "textbook.xml.gz").resolve())
        return read_sbml_model(textbook_sbml)
    elif model_name == "mini":
        mini_sbml = str((data_dir / "mini_fbc2.xml").resolve())
        return read_sbml_model(mini_sbml)
    elif model_name == "salmonella":
        salmonella_pickle = str((data_dir / "salmonella.pickle").resolve())
        model_name = salmonella_pickle
    with open(model_name, mode="rb") as infile:
        return _load(infile)


@pytest.fixture(scope="session")
def data_directory() -> Path:
    """Provide session-level fixture for test data directory."""
    return data_dir


@pytest.fixture(scope="session")
def cobra_data_directory() -> Path:
    """Provide session-level fixture for cobra data directory."""
    return Path(cobra_data_dir)


@pytest.fixture(scope="session")
def empty_once() -> Model:
    """Provide session-level fixture for empty model."""
    return Model()


@pytest.fixture(scope="function")
def empty_model(empty_once: Model) -> Model:
    """Provide function-level fixture for empty model."""
    return empty_once.copy()


@pytest.fixture(scope="session")
def small_model() -> Model:
    """Provide session-level fixture for textbook model."""
    return create_test_model("textbook")


@pytest.fixture(scope="function")
def model(small_model: Model) -> Model:
    """Provide function-level fixture for textbook model."""
    return small_model.copy()


@pytest.fixture(scope="session")
def large_once() -> Model:
    """Provide session-level fixture for ecoli model."""
    return create_test_model("ecoli")


@pytest.fixture(scope="function")
def large_model(large_once: Model) -> Model:
    """Provide function-level fixture for ecoli model."""
    return large_once.copy()


@pytest.fixture(scope="session")
def medium_model() -> Model:
    """Provide session-level fixture for salmonella model."""
    return create_test_model("salmonella")


@pytest.fixture(scope="function")
def salmonella(medium_model: Model) -> Model:
    """Provide function-level fixture for salmonella model."""
    return medium_model.copy()


@pytest.fixture(scope="function")
def solved_model(data_directory: Path) -> Tuple[Solution, Model]:
    """Provide function-level fixture for solved textbook model."""
    model = create_test_model("textbook")
    with (data_directory / "textbook_solution.pickle").open(mode="rb") as infile:
        solution = _load(infile)
    return solution, model


@pytest.fixture(scope="session")
def tiny_toy_model() -> Model:
    """Provide session-level fixture for tiny toy model."""
    tiny = Model("Toy Model")
    m1 = Metabolite("M1")
    d1 = Reaction("ex1")
    d1.add_metabolites({m1: -1})
    d1.upper_bound = 0
    d1.lower_bound = -1000
    tiny.add_reactions([d1])
    tiny.objective = "ex1"
    return tiny


stable_optlang = ["glpk", "cplex", "gurobi"]
all_solvers = ["optlang-" + s for s in stable_optlang if s in sutil.solvers]


@pytest.fixture(params=all_solvers, scope="session")
def opt_solver(request: pytest.FixtureRequest) -> str:
    """Provide session-level fixture for parametrized optlang solver names."""
    return request.param


@pytest.fixture(scope="function")
def metabolites(model: Model, request: pytest.FixtureRequest) -> List[Metabolite]:
    """Provide function-level fixture for metabolite set based on `request`."""
    if request.param == "exchange":
        return [
            met
            for met in model.metabolites
            if met.compartment == "e" and "EX_" + met.id not in model.reactions
        ]
    elif request.param == "demand":
        return [
            met
            for met in model.metabolites
            if met.compartment == "c" and "DM_" + met.id not in model.reactions
        ]
    elif request.param == "sink":
        return [
            met
            for met in model.metabolites
            if met.compartment == "c" and "SK_" + met.id not in model.reactions
        ]
    else:
        raise ValueError("Unknown metabolites {}".format(request.param))
