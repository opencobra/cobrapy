from __future__ import with_statement, absolute_import
import sys
from os import name as __name
available_tests = ['unit_tests', 'solvers', 'flux_analysis', 'io_tests']

del __name

from os.path import abspath as __abspath
from os.path import join as __join
from os.path import split as __split
from os.path import sep as __sep

cobra_directory = __abspath(__join(__split(__abspath(__file__))[0], ".."))
cobra_location = __abspath(__join(cobra_directory, ".."))
data_directory = __join(__split(__abspath(__file__))[0], "data")
if not data_directory.endswith(__sep):
    data_directory += __sep

salmonella_sbml = __join(data_directory, "salmonella.xml")
salmonella_fbc_sbml = __join(data_directory, "salmonella_fbc.xml")
salmonella_pickle = __join(data_directory, "salmonella.pickle")
salmonella_reaction_p_values_pickle = __join(data_directory, "salmonella_reaction_p_values.pickle")
ecoli_sbml = __join(data_directory, "iJO1366.xml")
ecoli_pickle = __join(data_directory, "iJO1366.pickle")
ecoli_mat = __join(data_directory, "iJO1366.mat")
ecoli_json = __join(data_directory, "iJO1366.json")
yersinia_sbml = __join(data_directory, 'Yersinia_pestis_CO92_iPC815.xml')
yersinia_pickle = __join(data_directory, 'Yersinia_pestis_CO92_iPC815.pickle')

__test_pickles = {'Salmonella_enterica': salmonella_pickle,
                  'Escherichia_coli': ecoli_pickle,
                  'Yersinia_pestis': yersinia_pickle}
__test_xml = {'Salmonella_enterica': salmonella_sbml,
              'Escherichia_coli': ecoli_sbml,
              'Yersinia_pestis': yersinia_sbml}
del __abspath, __join, __split, __sep

def create_test_model(test_pickle=salmonella_pickle):
    """Returns a cobra model for testing.  The default model is the up to date
    version of the Salmonella enterica Typhimurium LT2 model published in
    Thiele et al. 2011 BMC Sys Bio 5:8

    test_pickle: The complete file name of a pickled cobra.Model or SBML XML
    file to be read.  We currently provide Salmonella enterica Typhimurium
    and Escherichia coli models whose paths are stored in cobra.test.salmonella_pickle
    and cobra.test.ecoli_pickle, respectively.  The ecoli model is a variant of the
    model published in Orth et al. 2011 Mol Syst Biol 7:535

    """
    from os import name as __name
    try:
        from cPickle import load
    except:
        from pickle import load

    try: 
        with open(test_pickle, "rb") as infile:
            model = load(infile)
    except:
        #if the pickle can't be loaded then load the sbml xml
        from warnings import warn
        warn("Couldn't load %s.  Loading the default model %s instead"%(test_pickle,
                                                                        salmonella_sbml))
        sys.path.insert(0, cobra_location)
        from cobra.io import read_sbml_model
        model = read_sbml_model(salmonella_sbml)
        sys.path.pop(0)
    return model

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
