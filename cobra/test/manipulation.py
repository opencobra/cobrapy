from unittest import TestCase, TestLoader, TextTestRunner

import sys

if __name__ == "__main__":
    sys.path.insert(0, "../..")
    from cobra.test import create_test_model, data_directory
    from cobra.core import Metabolite
    from cobra.manipulation import *
    sys.path.pop(0)
else:
    from . import create_test_model, data_directory
    from ..core import Metabolite
    from ..manipulation import *


class TestManipulation(TestCase):
    """Test functions in cobra.manipulation"""

    def test_canonical_form(self):
        model = create_test_model("textbook")
        # add G constraint to test
        g_constr = Metabolite("SUCCt2_2__test_G_constraint")
        g_constr._constraint_sense = "G"
        g_constr._bound = 5.0
        model.reactions.get_by_id("SUCCt2_2").add_metabolites({g_constr: 1})
        self.assertAlmostEqual(model.optimize("maximize").f, 0.855, places=3)
        # convert to canonical form
        model = canonical_form(model)
        self.assertAlmostEqual(model.optimize("maximize").f, 0.855, places=3)

    def test_canonical_form_minimize(self):
        model = create_test_model("textbook")
        # make a minimization problem
        model.reactions.get_by_id("Biomass_Ecoli_core").lower_bound = 0.5
        for reaction in model.reactions:
            reaction.objective_coefficient = reaction.id == "GAPD"
        self.assertAlmostEqual(model.optimize("minimize").f, 6.27, places=3)
        # convert to canonical form. Convert minimize to maximize
        model = canonical_form(model, objective_sense="minimize")
        self.assertAlmostEqual(model.optimize("maximize").f, -6.27, places=3)
        # lower bounds should now be <= constraints
        self.assertEqual(
            model.reactions.get_by_id("Biomass_Ecoli_core").lower_bound, 0.0)


# make a test suite to run all of the tests
loader = TestLoader()
suite = loader.loadTestsFromModule(sys.modules[__name__])


def test_all():
    TextTestRunner(verbosity=2).run(suite)

if __name__ == "__main__":
    test_all()
