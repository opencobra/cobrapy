from __future__ import absolute_import
from os.path import join, abspath, dirname
import unittest as _unittest

try:
    from cPickle import load as _load
except:
    from pickle import load as _load

from ..io import read_sbml_model


available_tests = ['unit_tests', 'solvers', 'flux_analysis', 'io_tests',
                   'design', 'manipulation']


cobra_directory = abspath(join(dirname(abspath(__file__)), ".."))
cobra_location = abspath(join(cobra_directory, ".."))
data_directory = join(cobra_directory, "test", "data", "")

salmonella_sbml = join(data_directory, "salmonella.xml")
salmonella_pickle = join(data_directory, "salmonella.pickle")

ecoli_sbml = join(data_directory, "iJO1366.xml")
textbook_sbml = join(data_directory, "textbook.xml.gz")
mini_sbml = join(data_directory, "mini_fbc2.xml")

del abspath, join, dirname


def create_test_model(model_name="salmonella"):
    """Returns a cobra model for testing

    model_name: str
        One of 'ecoli', 'textbook', or 'salmonella', or the
        path to a pickled cobra.Model

    """

    if model_name == "ecoli":
        return read_sbml_model(ecoli_sbml)
    elif model_name == "textbook":
        return read_sbml_model(textbook_sbml)
    elif model_name == "mini":
        return read_sbml_model(mini_sbml)

    if model_name == "salmonella":
        model_name = salmonella_pickle
    with open(model_name, "rb") as infile:
        return _load(infile)


def create_test_suite():
    """create a unittest.TestSuite with available tests"""
    loader = _unittest.TestLoader()
    suite = _unittest.TestSuite()
    for test_name in available_tests:
        exec("from . import " + test_name)
        suite.addTests(loader.loadTestsFromModule(eval(test_name)))
    return suite

suite = create_test_suite()


def test_all():
    """###running unit tests on cobra py###"""
    status = not _unittest.TextTestRunner(verbosity=2).run(
        create_test_suite()
    ).wasSuccessful()
    return status
