from os.path import join, abspath, dirname
from cobra.io import read_sbml_model
import pytest

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


def test_all():
    """ alias for running all unit-tests on installed cobra
    """
    return pytest.main(
        ['--pyargs', 'cobra', '--benchmark-skip', '-v', '-rs']) == 0
