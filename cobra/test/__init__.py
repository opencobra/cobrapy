# -*- coding: utf-8 -*-

from __future__ import absolute_import

from os.path import abspath, dirname, join

from cobra.io import read_sbml_model

try:
    import pytest
    import pytest_benchmark
except ImportError:
    pytest = None
try:
    from cPickle import load as _load
except ImportError:
    from pickle import load as _load

cobra_directory = abspath(join(dirname(abspath(__file__)), ".."))
cobra_location = abspath(join(cobra_directory, ".."))
data_dir = join(cobra_directory, "test", "data", "")


def create_test_model(model_name="salmonella"):
    """Returns a cobra model for testing

    model_name: str
        One of 'ecoli', 'textbook', or 'salmonella', or the
        path to a pickled cobra.Model

    """
    if model_name == "ecoli":
        ecoli_sbml = join(data_dir, "iJO1366.xml")
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
    """ alias for running all unit-tests on installed cobra
    """
    if pytest:
        args = args if args else []

        return pytest.main(
            ['--pyargs', 'cobra', '--benchmark-skip', '-v', '-rs'] + args
        )
    else:
        raise ImportError('missing package pytest and pytest_benchmark'
                          ' required for testing')
