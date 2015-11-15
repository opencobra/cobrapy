from unittest import TestCase, TestLoader, TextTestRunner, skipIf

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
        self.assertAlmostEqual(model.optimize("Maximize").f, 0.874, places=3)
        dual = dual_problem(model)
        self.assertAlmostEqual(dual.optimize("Minimize").f, 0.874, places=3)

    def test_dual_integer_vars_as_lp(self):
        model = create_test_model("textbook")
        var = _add_decision_variable(model, "GAPD")
        self.assertAlmostEqual(model.optimize("Maximize").f, 0.874, places=3)
        # as lp: make integer continuous, set to 1
        dual = dual_problem(model, "Maximize", [var.id], copy=True)
        r = dual.reactions.get_by_id(var.id)
        r.lower_bound = r.upper_bound = 1
        r.variable_kind = "continuous"
        self.assertAlmostEqual(dual.optimize("Minimize").f, 0.874, places=3)

    @skipIf(no_mip_solver, "no MILP solver found")
    def test_dual_integer_vars_as_mip(self):
        try:
            solver = get_solver_name(mip=True)
        except:
            raise Exception("no MILP solver found")

        model = create_test_model("textbook")
        var = _add_decision_variable(model, "GAPD")
        dual = dual_problem(model, "Maximize", [var.id], copy=True)
        var_in_dual = dual.reactions.get_by_id(var.id)

        # turn off GAPD in model
        var.lower_bound = var.upper_bound = 0
        self.assertAlmostEqual(model.optimize("Maximize", presolve=True).f, 0.0, places=3)
        # should not affect copied dual problem
        self.assertAlmostEqual(dual.optimize("Minimize", presolve=True).f, 0.874, places=3)

        # turn off GAPD in dual
        var_in_dual.lower_bound = var_in_dual.upper_bound = 1
        self.assertAlmostEqual(dual.optimize("Minimize", presolve=True).f, 0.874, places=3)

        # turn on GAPD in dual
        var_in_dual.lower_bound = var_in_dual.upper_bound = 0
        self.assertAlmostEqual(dual.optimize("Minimize", presolve=True).f, 0.0, places=3)

    @skipIf(no_mip_solver, "no MILP solver found")
    def test_optknock(self):
        model = create_test_model("textbook")
        knockable_reactions = ["ACKr", "AKGDH", "ACALD", "LDH_L"]
        optknock_problem = set_up_optknock(model, "EX_lac__L_e",
                                           model.reactions, n_knockouts=2,
                                           copy=False)
        solution = run_optknock()
        self.assertIn("ACKr", solution.knockouts)
        self.assertIn("ACALD", solution.knockouts)

# make a test suite to run all of the tests
loader = TestLoader()
suite = loader.loadTestsFromModule(sys.modules[__name__])


def test_all():
    TextTestRunner(verbosity=2).run(suite)

if __name__ == "__main__":
    test_all()
