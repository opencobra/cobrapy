from __future__ import with_statement, absolute_import
import sys
available_tests = ['unit_tests', 'solvers', 'flux_analysis', 'io_tests']

from os.path import abspath, join, split, sep

cobra_directory = abspath(join(split(abspath(__file__))[0], ".."))
cobra_location = abspath(join(cobra_directory, ".."))
data_directory = join(split(abspath(__file__))[0], "data")
if not data_directory.endswith(sep):
    data_directory += sep

salmonella_sbml = join(data_directory, "salmonella.xml")
salmonella_fbc_sbml = join(data_directory, "salmonella_fbc.xml")
salmonella_pickle = join(data_directory, "salmonella.pickle")
ecoli_sbml = join(data_directory, "iJO1366.xml")
ecoli_pickle = join(data_directory, "iJO1366.pickle")
ecoli_mat = join(data_directory, "iJO1366.mat")
ecoli_json = join(data_directory, "iJO1366.json")

__test_pickles = {'Salmonella_enterica': salmonella_pickle,
                  'Escherichia_coli': ecoli_pickle,
                  }
__test_xml = {'Salmonella_enterica': salmonella_sbml,
              'Escherichia_coli': ecoli_sbml,
              }
del abspath, join, split, sep


def create_test_model(test_pickle=salmonella_pickle):
    """Returns a cobra model for testing.  The default model is the up to date
    version of the Salmonella enterica Typhimurium LT2 model published in
    Thiele et al. 2011 BMC Sys Bio 5:8

    test_pickle: The filename of a pickled cobra.Model
    We currently provide Salmonella enterica Typhimurium and Escherichia coli
    models whose paths are stored in cobra.test.salmonella_pickle and
    cobra.test.ecoli_pickle, respectively.  The ecoli model is a variant of the
    model published in Orth et al. 2011 Mol Syst Biol 7:535

    """
    try:
        from cPickle import load
    except:
        from pickle import load

    if test_pickle == "salmonella":
        test_pickle = salmonella_pickle
    elif test_pickle == "ecoli":
        test_pickle = ecoli_pickle
    with open(test_pickle, "rb") as infile:
        return load(infile)


def create_test_suite():
    """create a unittest.TestSuite with available tests"""
    from unittest import TestLoader, TestSuite
    loader = TestLoader()
    suite = TestSuite()
    for test_name in available_tests:
        exec("from . import " + test_name)
        suite.addTests(loader.loadTestsFromModule(eval(test_name)))
    return suite

suite = create_test_suite()


def test_all():
    """###running unit tests on cobra py###"""
    from unittest import TextTestRunner
    TextTestRunner(verbosity=2).run(create_test_suite())
