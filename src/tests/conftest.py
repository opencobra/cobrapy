"""Define global fixtures."""

from os.path import abspath, dirname, join
from pickle import load as _load

from cobra import Metabolite, Model, Reaction
from cobra.io import read_sbml_model
from cobra.util import solver as sutil


try:
    import pytest
    import pytest_benchmark
except ImportError:
    pytest = None


cobra_directory = abspath(join(dirname(abspath(__file__)), ".."))
cobra_location = abspath(join(cobra_directory, ".."))
data_dir = join(cobra_directory, "tests", "data", "")


def create_test_model(model_name="salmonella"):
    """Return a cobra model for testing.

    model_name: str
        One of 'ecoli', 'textbook', or 'salmonella', or the
        path to a pickled cobra.Model

    """
    if model_name == "ecoli":
        ecoli_sbml = join(data_dir, "iJO1366.xml.gz")
        return read_sbml_model(ecoli_sbml)
    elif model_name == "textbook":
        textbook_sbml = join(data_dir, "textbook.xml.gz")
        return read_sbml_model(textbook_sbml)
    elif model_name == "mini":
        mini_sbml = join(data_dir, "mini_fbc2.xml")
        return read_sbml_model(mini_sbml)
    elif model_name == "salmonella":
        salmonella_pickle = join(data_dir, "salmonella.pickle")
        model_name = salmonella_pickle
    with open(model_name, "rb") as infile:
        return _load(infile)


def test_all(args=None):
    """Alias for running all unit-tests on installed cobra."""
    if pytest:
        args = args if args else []

        return pytest.main(
            ["--pyargs", "cobra", "--benchmark-skip", "-v", "-rs"] + args
        )
    else:
        raise ImportError(
            "missing package pytest and pytest_benchmark required for testing"
        )


def pytest_addoption(parser):
    try:
        parser.addoption("--run-slow", action="store_true", help="run slow tests")
        parser.addoption(
            "--run-non-deterministic",
            action="store_true",
            help="run tests that sometimes (rarely) fail",
        )
    except ValueError:
        pass


@pytest.fixture(scope="session")
def data_directory():
    return data_dir


@pytest.fixture(scope="session")
def empty_once():
    return Model()


@pytest.fixture(scope="function")
def empty_model(empty_once):
    return empty_once.copy()


@pytest.fixture(scope="session")
def small_model():
    return create_test_model("textbook")


@pytest.fixture(scope="function")
def model(small_model):
    return small_model.copy()


@pytest.fixture(scope="session")
def large_once():
    return create_test_model("ecoli")


@pytest.fixture(scope="function")
def large_model(large_once):
    return large_once.copy()


@pytest.fixture(scope="session")
def medium_model():
    return create_test_model("salmonella")


@pytest.fixture(scope="function")
def salmonella(medium_model):
    return medium_model.copy()


@pytest.fixture(scope="function")
def solved_model(data_directory):
    model = create_test_model("textbook")
    with open(join(data_directory, "textbook_solution.pickle"), "rb") as infile:
        solution = _load(infile)
    return solution, model


@pytest.fixture(scope="session")
def tiny_toy_model():
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
def opt_solver(request):
    return request.param


@pytest.fixture(scope="function")
def metabolites(model, request):
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
        raise ValueError("unknown metabolites {}".format(request.param))
