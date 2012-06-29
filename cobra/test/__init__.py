from os import name as __name
from sys import modules as __modules
from warnings import warn
if __name == 'java':
    raise Exception("%s is not yet supported on jython"%__modules[__name__])

del __name, __modules

from os.path import abspath as __abspath
from os.path import join as __join
from os.path import split as __split
from os.path import sep as __sep

data_directory = __join(__split(__abspath(__file__))[0], "data")
if not data_directory.endswith(__sep):
    data_directory += __sep

salmonella_sbml = __join(data_directory, "salmonella.xml")
salmonella_pickle = __join(data_directory, "salmonella.pickle")
salmonella_reaction_p_values_pickle = __join(data_directory, "salmonella_reaction_p_values.pickle")
ecoli_sbml = __join(data_directory, "iJO1366.xml")
ecoli_pickle = __join(data_directory, "iJO1366.pickle")
ecoli_mat = __join(data_directory, "iJO1366.mat")

def create_test_model(test_pickle=salmonella_pickle):
    """returns the salmonella model for testing"""
    try:
        from cPickle import load
    except:
        from pickle import load
    with open(test_pickle, "rb") as infile:
        model = load(infile)
    return model


def test_all():
    """###running unit tests on cobra py###"""
    import sys
    sys.path.insert(0, "../..")
    import unittest
    from cobra.test.unit_tests import test_all
    print '###running general unit tests###'
    test_all()
    from cobra.test.flux_analysis import test_all
    print '\n###running flux_analysis unit tests###'
    test_all()
    print '\n###running solver unit tests###'
    from cobra.test.solvers import test_all
    test_all()
    sys.path.pop(0)

del __abspath, __join, __split, __sep

if __name__ == '__main__':
    test_all()
