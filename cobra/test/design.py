from unittest import TestCase, TestLoader, TextTestRunner

import sys

if __name__ == "__main__":
    sys.path.insert(0, "../..")
    from cobra.test import create_test_model, data_directory
    from cobra.design import *
    from cobra.design.design_algorithms import _add_decision_variable
    sys.path.pop(0)
else:
    from . import create_test_model, data_directory
    from ..design import *
    from ..design.design_algorithms import _add_decision_variable


class TestDesignAlgorithms(TestCase):
    """Test functions in cobra.design"""

    def test_dual(self):
        model = create_test_model("textbook")
        self.assertAlmostEqual(model.optimize("Maximize").f, 0.874, places=3)
        dual = dual_problem(model)
        self.assertAlmostEqual(dual.optimize("Minimize").f, 0.874, places=3)

    def test_dual_integer_vars(self):
        model = create_test_model("textbook")
        var = _add_decision_variable(model, "GAPD")
        self.assertAlmostEqual(model.optimize("Maximize").f, 0.874, places=3)
        dual = dual_problem(model, [var.id])
        self.assertAlmostEqual(dual.optimize("Minimize").f, 0.874, places=3)

    def test_dual_integer_vars_knockout(self):
        model = create_test_model("textbook")
        var = _add_decision_variable(model, "GAPD")
        dual = dual_problem(model, [var.id], copy=True)
        # turn off GAPD in model
        var.upper_bound = 0
        self.assertAlmostEqual(model.optimize("Maximize").f, 0.0, places=3)
        self.assertAlmostEqual(dual.optimize("Minimize").f, 0.874, places=3)
        # turn off GAPD in dual
        dual.reactions.get_by_id("GAPD_decision_var").upper_bound = 0
        self.assertAlmostEqual(dual.optimize("Minimize").f, 0.0, places=3)



# make a test suite to run all of the tests
loader = TestLoader()
suite = loader.loadTestsFromModule(sys.modules[__name__])


def test_all():
    TextTestRunner(verbosity=2).run(suite)

if __name__ == "__main__":
    test_all()
