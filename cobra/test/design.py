from unittest import TestCase, TestLoader, TextTestRunner

import sys

if __name__ == "__main__":
    sys.path.insert(0, "../..")
    from cobra.test import create_test_model, data_directory
    from cobra.design import *
    sys.path.pop(0)
else:
    from . import create_test_model, data_directory
    from ..design import *


class TestDesignAlgorithms(TestCase):
    """Test functions in cobra.design"""

    def test_dual(self):
        model = create_test_model("textbook")
        dual = dual_problem(model)
        self.assertAlmostEqual(dual.optimize('Minimize').f,
                               model.optimize('Maximize').f,
                               places=3)


# make a test suite to run all of the tests
loader = TestLoader()
suite = loader.loadTestsFromModule(sys.modules[__name__])


def test_all():
    TextTestRunner(verbosity=2).run(suite)

if __name__ == "__main__":
    test_all()
