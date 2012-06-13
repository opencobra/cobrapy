from os.path import abspath as _abspath
from os.path import join as _join
from os.path import split as _split
from os.path import sep as _sep

data_directory = _join(_split(_abspath(__file__))[0], "data")
if not data_directory.endswith(_sep):
    data_directory += _sep

salmonella_sbml = _join(data_directory, "salmonella.xml")
salmonella_pickle = _join(data_directory, "salmonella.pickle")
salmonella_reaction_p_values_pickle = _join(data_directory, "salmonella_reaction_p_values.pickle")
ecoli_sbml = _join(data_directory, "iJO1366.xml")
ecoli_pickle = _join(data_directory, "iJO1366.pickle")
ecoli_mat = _join(data_directory, "iJO1366.mat")

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
    """run unittests on the cobra toolbox"""
    import unittest
    from cobra.test.unit_tests import test_all
    print 'running general unit_tests'
    test_all()

    #print 'running new solver unit tests'
    #from cobra.test.solvers import test_all
    #test_all()


del _abspath
del _join
del _split
del _sep
