from unittest import TestCase, TestLoader, TextTestRunner, skipIf

import sys

if __name__ == "__main__":
    sys.path.insert(0, "../..")
    from cobra.test import create_test_model, data_directory
    from cobra.design import *
    from cobra.design.design_algorithms import _add_decision_variable
    from cobra.solvers import get_solver_name
    sys.path.pop(0)
else:
    from . import create_test_model, data_directory
    from ..design import *
    from ..design.design_algorithms import _add_decision_variable
    from ..solvers import get_solver_name

try:
    solver = get_solver_name(mip=True)
except:
    no_mip_solver = True
else:
    no_mip_solver = False


class TestDesignAlgorithms(TestCase):
    """Test functions in cobra.design"""

    def test_dual(self):
        model = create_test_model("textbook")
        self.assertAlmostEqual(model.optimize("maximize").f, 0.874, places=3)
        dual = dual_problem(model)
        self.assertAlmostEqual(dual.optimize("minimize").f, 0.874, places=3)

    def test_dual_integer_vars_as_lp(self):
        model = create_test_model("textbook")
        var = _add_decision_variable(model, "AKGDH")
        self.assertAlmostEqual(model.optimize("maximize").f, 0.874, places=3)
        # as lp: make integer continuous, set to 1
        dual = dual_problem(model, "maximize", [var.id], copy=True)
        r = dual.reactions.get_by_id(var.id)
        r.variable_kind = "continuous"
        r.lower_bound = r.upper_bound = 1
        self.assertAlmostEqual(dual.optimize("minimize").f, 0.874, places=3)
        r.lower_bound = r.upper_bound = 0
        self.assertAlmostEqual(dual.optimize("minimize").f, 0.858, places=3)

    @skipIf(no_mip_solver, "no MILP solver found")
    def test_dual_integer_vars_as_mip(self):
        # mip
        model = create_test_model("textbook")
        var = _add_decision_variable(model, "AKGDH")
        dual = dual_problem(model, "maximize", [var.id], copy=True)
        var_in_dual = dual.reactions.get_by_id(var.id)

        # minimization, so the optimal value state is to turn off AKGDH
        self.assertAlmostEqual(dual.optimize("minimize").f, 0.858, places=3)

        # turn off AKGDH in dual
        var_in_dual.lower_bound = var_in_dual.upper_bound = 1
        self.assertAlmostEqual(dual.optimize("minimize").f, 0.874, places=3)

        # turn on AKGDH in dual
        var_in_dual.lower_bound = var_in_dual.upper_bound = 0
        self.assertAlmostEqual(dual.optimize("minimize").f, 0.858, places=3)

    @skipIf(no_mip_solver, "no MILP solver found")
    def test_optknock(self):
        model = create_test_model("textbook")
        model.reactions.get_by_id("EX_o2_e").lower_bound = 0
        knockable_reactions = ["ACKr", "AKGDH", "ACALD", "LDH_D"]
        optknock_problem = set_up_optknock(model, "EX_lac__D_e",
                                           knockable_reactions, n_knockouts=2,
                                           copy=False)
        solution = run_optknock(optknock_problem, tolerance_integer=1e-9)
        self.assertIn("ACKr", solution.knockouts)
        self.assertIn("ACALD", solution.knockouts)
        self.assertAlmostEqual(solution.f, 17.891, places=3)

# make a test suite to run all of the tests
loader = TestLoader()
suite = loader.loadTestsFromModule(sys.modules[__name__])


def test_all():
    TextTestRunner(verbosity=2).run(suite)

if __name__ == "__main__":
    test_all()
