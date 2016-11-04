from os.path import join
from . import create_test_model, data_dir
import pytest
try:
    from cPickle import load as _load
except ImportError:
    from pickle import load as _load
import json


@pytest.fixture(scope="session")
def data_directory():
    return data_dir


@pytest.fixture(scope="function")
def model():
    return create_test_model("textbook")


@pytest.fixture(scope="function")
def large_model():
    return create_test_model("ecoli")


@pytest.fixture(scope="function")
def array_model():
    return create_test_model("textbook").to_array_based_model()


@pytest.fixture(scope="function")
def salmonella():
    return create_test_model("salmonella")


@pytest.fixture(scope="function")
def solved_model(data_directory):
    model = create_test_model("textbook")
    with open(join(data_directory, "textbook_solution.pickle"),
              "rb") as infile:
        model.solution = _load(infile)
    return model


@pytest.fixture(scope="function")
def fva_results(data_directory):
    with open(join(data_directory, "textbook_fva.json"), "r") as infile:
        return json.load(infile)
