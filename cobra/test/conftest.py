# -*- coding: utf-8 -*-
from __future__ import absolute_import

import json
from os.path import join

import pytest

from . import create_test_model, data_dir

try:
    from cPickle import load as _load
except ImportError:
    from pickle import load as _load


def pytest_addoption(parser):
    try:
        parser.addoption("--run-slow", action="store_true",
                         help="run slow tests")
        parser.addoption("--run-non-deterministic", action="store_true",
                         help="run tests that sometimes (rarely) fail")
    except ValueError:
        pass


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
