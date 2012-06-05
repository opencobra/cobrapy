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


def run_tests():
    """run unittests on the cobra toolbox"""
    import unittest
    # NOTE: This only runs a test suite in the oven
    from cobra.test.unit_tests import suite as test_suite
    # if/when another test suite is written, it can also get tested
    unittest.TextTestRunner(verbosity=2).run(test_suite)
    

if __name__ == "__main__":
    run_tests()
