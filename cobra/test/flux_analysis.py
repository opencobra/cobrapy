from unittest import TestCase, TestLoader, TextTestRunner, skipIf

from warnings import warn
import sys
from os.path import join
from os import name
from json import load
from contextlib import contextmanager
import pickle
import re

from six import iteritems, StringIO

try:
    import numpy
except:
    numpy = None
try:
    import matplotlib
except:
    matplotlib = None
try:
    import pandas
except:
    pandas = None
try:
    import tabulate
except:
    tabulate = None

if __name__ == "__main__":
    sys.path.insert(0, "../..")
    from cobra.test import create_test_model, data_directory
    from cobra import Model, Reaction, Metabolite
    from cobra.manipulation import initialize_growth_medium
    from cobra.solvers import solver_dict, get_solver_name
    from cobra.flux_analysis import *
    sys.path.pop(0)
else:
    from . import create_test_model, data_directory
    from .. import Model, Reaction, Metabolite
    from ..manipulation import initialize_growth_medium
    from ..solvers import solver_dict, get_solver_name
    from ..flux_analysis import *


@contextmanager
def captured_output():
    """ A context manager to test the IO summary methods """
    new_out, new_err = StringIO(), StringIO()
    old_out, old_err = sys.stdout, sys.stderr
    try:
        sys.stdout, sys.stderr = new_out, new_err
        yield sys.stdout, sys.stderr
    finally:
        sys.stdout, sys.stderr = old_out, old_err


class TestCobraFluxAnalysis(TestCase):
    """Test the simulation functions in cobra.flux_analysis"""

    def setUp(self):
        pass

    def test_pFBA(self):
        model = create_test_model("textbook")
        for solver in solver_dict:
            optimize_minimal_flux(model, solver=solver)
            abs_x = [abs(i) for i in model.solution.x]
            self.assertEqual(model.solution.status, "optimal")
            self.assertAlmostEqual(model.solution.f, 0.8739, places=3)
            self.assertAlmostEqual(sum(abs_x), 518.4221, places=3)

            # Test desired_objective_value
            desired_objective = 0.8
            optimize_minimal_flux(model, solver=solver,
                                  desired_objective_value=desired_objective)
            abs_x = [abs(i) for i in model.solution.x]
            self.assertEqual(model.solution.status, "optimal")
            self.assertAlmostEqual(model.solution.f, desired_objective,
                                   places=3)
            self.assertAlmostEqual(sum(abs_x), 476.1594, places=3)

            # Test fraction_of_optimum
            optimize_minimal_flux(model, solver=solver,
                                  fraction_of_optimum=0.95)
            abs_x = [abs(i) for i in model.solution.x]
            self.assertEqual(model.solution.status, "optimal")
            self.assertAlmostEqual(model.solution.f, 0.95*0.8739, places=3)
            self.assertAlmostEqual(sum(abs_x), 493.4400, places=3)

            # Make sure the model works for non-unity objective values
            model.reactions.Biomass_Ecoli_core.objective_coefficient = 2
            optimize_minimal_flux(model, solver=solver)
            self.assertAlmostEqual(model.solution.f, 2*0.8739, places=3)
            model.reactions.Biomass_Ecoli_core.objective_coefficient = 1

            # Test some erroneous inputs -- multiple objectives
            model.reactions.ATPM.objective_coefficient = 1
            with self.assertRaises(ValueError):
                optimize_minimal_flux(model, solver=solver)
            model.reactions.ATPM.objective_coefficient = 0

            # Minimization of objective
            with self.assertRaises(ValueError):
                optimize_minimal_flux(model, solver=solver,
                                      objective_sense='minimize')

            # Infeasible solution
            atpm = float(model.reactions.ATPM.lower_bound)
            model.reactions.ATPM.lower_bound = 500
            with self.assertRaises(ValueError):
                optimize_minimal_flux(model, solver=solver)
            model.reactions.ATPM.lower_bound = atpm

    def test_single_gene_deletion_fba(self):
        cobra_model = create_test_model("textbook")
        # expected knockouts for textbook model
        growth_dict = {"b0008": 0.87, "b0114": 0.80, "b0116": 0.78,
                       "b2276": 0.21, "b1779": 0.00}

        rates, statuses = single_gene_deletion(cobra_model,
                                               gene_list=growth_dict.keys(),
                                               method="fba")
        for gene, expected_value in iteritems(growth_dict):
            self.assertEqual(statuses[gene], 'optimal')
            self.assertAlmostEqual(rates[gene], expected_value, places=2)

    def test_single_gene_deletion_moma(self):
        # MOMA requires a QP solver
        try:
            get_solver_name(qp=True)
        except:
            self.skipTest("no qp support")

        cobra_model = create_test_model("textbook")
        # expected knockouts for textbook model
        growth_dict = {"b0008": 0.87, "b0114": 0.71, "b0116": 0.56,
                       "b2276": 0.11, "b1779": 0.00}

        rates, statuses = single_gene_deletion(cobra_model,
                                               gene_list=growth_dict.keys(),
                                               method="moma")
        for gene, expected_value in iteritems(growth_dict):
            self.assertEqual(statuses[gene], 'optimal')
            self.assertAlmostEqual(rates[gene], expected_value, places=2)

    def test_single_reaction_deletion(self):
        cobra_model = create_test_model("textbook")
        expected_results = {'FBA': 0.70404, 'FBP': 0.87392, 'CS': 0,
                            'FUM': 0.81430, 'GAPD': 0, 'GLUDy': 0.85139}

        results, status = single_reaction_deletion(
            cobra_model, reaction_list=expected_results.keys())
        self.assertEqual(len(results), 6)
        self.assertEqual(len(status), 6)
        for status_value in status.values():
            self.assertEqual(status_value, "optimal")
        for reaction, value in results.items():
            self.assertAlmostEqual(value, expected_results[reaction], 5)

    def compare_matrices(self, matrix1, matrix2, places=3):
        nrows = len(matrix1)
        ncols = len(matrix1[0])
        self.assertEqual(nrows, len(matrix2))
        self.assertEqual(ncols, len(matrix2[0]))
        for i in range(nrows):
            for j in range(ncols):
                self.assertAlmostEqual(matrix1[i][j], matrix2[i][j],
                                       places=places)

    @skipIf(numpy is None, "double deletions require numpy")
    def test_double_gene_deletion(self):
        cobra_model = create_test_model("textbook")
        genes = ["b0726", "b4025", "b0724", "b0720",
                 "b2935", "b2935", "b1276", "b1241"]
        growth_list = [
            [0.858, 0.857, 0.814, 0.000, 0.858, 0.858, 0.858, 0.858],
            [0.857, 0.863, 0.739, 0.000, 0.863, 0.863, 0.863, 0.863],
            [0.814, 0.739, 0.814, 0.000, 0.814, 0.814, 0.814, 0.814],
            [0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000],
            [0.858, 0.863, 0.814, 0.000, 0.874, 0.874, 0.874, 0.874],
            [0.858, 0.863, 0.814, 0.000, 0.874, 0.874, 0.874, 0.874],
            [0.858, 0.863, 0.814, 0.000, 0.874, 0.874, 0.874, 0.874],
            [0.858, 0.863, 0.814, 0.000, 0.874, 0.874, 0.874, 0.874]]
        opts = {"number_of_processes": 1} if name == "nt" else {}
        solution = double_gene_deletion(cobra_model, gene_list1=genes, **opts)
        self.assertEqual(solution["x"], genes)
        self.assertEqual(solution["y"], genes)
        self.compare_matrices(growth_list, solution["data"])
        # test when lists differ slightly
        solution = double_gene_deletion(cobra_model, gene_list1=genes[:-1],
                                        gene_list2=genes,
                                        number_of_processes=1)
        self.assertEqual(solution["x"], genes[:-1])
        self.assertEqual(solution["y"], genes)
        self.compare_matrices(growth_list[:-1], solution["data"])

    @skipIf(numpy is None, "double deletions require numpy")
    def test_double_reaction_deletion(self):
        cobra_model = create_test_model("textbook")
        reactions = ['FBA', 'ATPS4r', 'ENO', 'FRUpts2']
        growth_list = [[0.704, 0.135, 0.000, 0.704],
                       [0.135, 0.374, 0.000, 0.374],
                       [0.000, 0.000, 0.000, 0.000],
                       [0.704, 0.374, 0.000, 0.874]]

        solution = double_reaction_deletion(cobra_model,
                                            reaction_list1=reactions,
                                            number_of_processes=1)
        self.assertEqual(solution["x"], reactions)
        self.assertEqual(solution["y"], reactions)
        self.compare_matrices(growth_list, solution["data"])

    def test_flux_variability(self):
        with open(join(data_directory, "textbook_fva.json"), "r") as infile:
            fva_results = load(infile)

        infeasible_model = create_test_model("textbook")
        infeasible_model.reactions.get_by_id("EX_glc__D_e").lower_bound = 0
        for solver in solver_dict:
            # esolver is really slow
            if solver == "esolver":
                continue
            cobra_model = create_test_model("textbook")
            fva_out = flux_variability_analysis(
                cobra_model, solver=solver,
                reaction_list=cobra_model.reactions[1::3])
            for name, result in iteritems(fva_out):
                for k, v in iteritems(result):
                    self.assertAlmostEqual(fva_results[name][k], v, places=5)

            # ensure that an infeasible model does not run FVA
            self.assertRaises(ValueError, flux_variability_analysis,
                              infeasible_model, solver=solver)

    def test_find_blocked_reactions(self):
        m = create_test_model("textbook")
        result = find_blocked_reactions(m, m.reactions[40:46])
        self.assertEqual(result, ['FRUpts2'])

        result = find_blocked_reactions(m, m.reactions[42:48])
        self.assertEqual(set(result), {'FUMt2_2', 'FRUpts2'})

        result = find_blocked_reactions(m, m.reactions[30:50],
                                        open_exchanges=True)
        self.assertEqual(result, [])

    def test_loopless(self):
        try:
            solver = get_solver_name(mip=True)
        except:
            self.skipTest("no MILP solver found")
        test_model = Model()
        test_model.add_metabolites(Metabolite("A"))
        test_model.add_metabolites(Metabolite("B"))
        test_model.add_metabolites(Metabolite("C"))
        EX_A = Reaction("EX_A")
        EX_A.add_metabolites({test_model.metabolites.A: 1})
        DM_C = Reaction("DM_C")
        DM_C.add_metabolites({test_model.metabolites.C: -1})
        v1 = Reaction("v1")
        v1.add_metabolites({test_model.metabolites.A: -1,
                            test_model.metabolites.B: 1})
        v2 = Reaction("v2")
        v2.add_metabolites({test_model.metabolites.B: -1,
                            test_model.metabolites.C: 1})
        v3 = Reaction("v3")
        v3.add_metabolites({test_model.metabolites.C: -1,
                            test_model.metabolites.A: 1})
        DM_C.objective_coefficient = 1
        test_model.add_reactions([EX_A, DM_C, v1, v2, v3])
        feasible_sol = construct_loopless_model(test_model).optimize()
        v3.lower_bound = 1
        infeasible_sol = construct_loopless_model(test_model).optimize()
        self.assertEqual(feasible_sol.status, "optimal")
        self.assertEqual(infeasible_sol.status, "infeasible")

    def test_gapfilling(self):
        try:
            solver = get_solver_name(mip=True)
        except:
            self.skipTest("no MILP solver found")
        m = Model()
        m.add_metabolites(map(Metabolite, ["a", "b", "c"]))
        r = Reaction("EX_A")
        m.add_reaction(r)
        r.add_metabolites({m.metabolites.a: 1})
        r = Reaction("r1")
        m.add_reaction(r)
        r.add_metabolites({m.metabolites.b: -1, m.metabolites.c: 1})
        r = Reaction("DM_C")
        m.add_reaction(r)
        r.add_metabolites({m.metabolites.c: -1})
        r.objective_coefficient = 1

        U = Model()
        r = Reaction("a2b")
        U.add_reaction(r)
        r.build_reaction_from_string("a --> b", verbose=False)
        r = Reaction("a2d")
        U.add_reaction(r)
        r.build_reaction_from_string("a --> d", verbose=False)

        # GrowMatch
        result = gapfilling.growMatch(m, U)[0]
        self.assertEqual(len(result), 1)
        self.assertEqual(result[0].id, "a2b")
        # SMILEY
        result = gapfilling.SMILEY(m, "b", U)[0]
        self.assertEqual(len(result), 1)
        self.assertEqual(result[0].id, "a2b")

        # 2 rounds of GrowMatch with exchange reactions
        result = gapfilling.growMatch(m, None, ex_rxns=True, iterations=2)
        self.assertEqual(len(result), 2)
        self.assertEqual(len(result[0]), 1)
        self.assertEqual(len(result[1]), 1)
        self.assertEqual({i[0].id for i in result},
                         {"SMILEY_EX_b", "SMILEY_EX_c"})

    @skipIf(numpy is None, "phase plane requires numpy")
    def test_phenotype_phase_plane(self):
        model = create_test_model("textbook")
        data = calculate_phenotype_phase_plane(
            model, "EX_glc__D_e", "EX_o2_e",
            reaction1_npoints=20, reaction2_npoints=20)
        self.assertEqual(data.growth_rates.shape, (20, 20))
        self.assertAlmostEqual(data.growth_rates.max(), 1.20898, places=4)
        self.assertAlmostEqual(abs(data.growth_rates[0, :]).max(), 0, places=4)
        if matplotlib is None:
            self.skipTest("can't test plots without matplotlib")
        data.plot()

    def check_entries(self, out, desired_entries):
        """ensure each entry in desired_entries appears in output"""
        output = out.getvalue().strip()
        output_set = set((re.sub('\s', '', l) for l in output.splitlines()))

        for item in desired_entries:
            self.assertIn(re.sub('\s', '', item), output_set)

    @skipIf((pandas is None) or (tabulate is None),
            "summary methods require pandas and tabulate")
    def test_summary_methods(self):

        # Test model summary methods
        model = create_test_model("textbook")
        with self.assertRaises(Exception):
            model.summary()

        # Load model solution
        with open(join(data_directory, "textbook_solution.pickle"),
                  "rb") as infile:
            model.solution = pickle.load(infile)

        desired_entries = [
            'idFluxRangeidFluxRangeBiomass_Ecol...0.874',
            'o2_e       21.8   [19.9, 23.7]'
            'h2o_e       29.2  [25, 30.7]',
            'glc__D_e   10     [9.52, 10]'
            'co2_e       22.8  [18.9, 24.7]',
            'nh4_e       4.77  [4.53, 5.16]'
            'h_e         17.5  [16.7, 22.4]',
            'pi_e        3.21  [3.05, 3.21]'
            'for_e        0    [0, 5.72]',
            'ac_e         0    [0, 1.91]',
            'pyr_e        0    [0, 1.27]',
            'lac__D_e     0    [0, 1.07]',
            'succ_e       0    [0, 0.837]',
            'glu__L_e     0    [0, 0.636]',
            'akg_e        0    [0, 0.715]',
            'etoh_e       0    [0, 1.11]',
            'acald_e      0    [0, 1.27]',
        ]
        for solver in solver_dict:
            with captured_output() as (out, err):
                model.summary(fva=0.95, solver=solver)
            self.check_entries(out, desired_entries)

        # test non-fva version (these should be fixed for textbook model
        desired_entries = [
            'o2_e      21.8',
            'glc__D_e  10',
            'nh4_e      4.77',
            'pi_e       3.21',
            'h2o_e  29.2',
            'co2_e  22.8',
            'h_e    17.5',
            'Biomass_Ecol...  0.874',
        ]
        # Need to use a different method here because
        # there are multiple entries per line.
        for solver in solver_dict:
            with captured_output() as (out, err):
                model.summary()

            s = out.getvalue()
            for i in desired_entries:
                self.assertIn(i, s)

        # Test metabolite summary methods
        desired_entries = [
            'PRODUCING REACTIONS -- Ubiquinone-8 (q8_c)',
            '%       FLUX  RXN ID    REACTION',
            '100%   43.6   CYTBD     '
            '2.0 h_c + 0.5 o2_c + q8h2_c --> h2o_c + 2.0 h_e...',
            'CONSUMING REACTIONS -- Ubiquinone-8 (q8_c)',
            '%       FLUX  RXN ID    REACTION',
            '88%    38.5   NADH16    '
            '4.0 h_c + nadh_c + q8_c --> 3.0 h_e + nad_c + q...',
            '12%     5.06  SUCDi     q8_c + succ_c --> fum_c + q8h2_c',
        ]

        for solver in solver_dict:
            with captured_output() as (out, err):
                model.metabolites.q8_c.summary()
            self.check_entries(out, desired_entries)

        desired_entries = [
            'PRODUCING REACTIONS -- D-Fructose 1,6-bisphosphate (fdp_c)',
            '----------------------------------------------------------',
            '%       FLUX  RANGE         RXN ID    REACTION',
            '100%    7.48  [6.17, 9.26]  PFK       '
            'atp_c + f6p_c --> adp_c + fdp_c + h_c',
            'CONSUMING REACTIONS -- D-Fructose 1,6-bisphosphate (fdp_c)',
            '----------------------------------------------------------',
            '%       FLUX  RANGE         RXN ID    REACTION',
            '100%    7.48  [6.17, 8.92]  FBA       fdp_c <=> dhap_c + g3p_c',
            '0%      0     [0, 1.72]     FBP       '
            'fdp_c + h2o_c --> f6p_c + pi_c',
        ]

        for solver in solver_dict:
            with captured_output() as (out, err):
                model.metabolites.fdp_c.summary(fva=0.99, solver=solver)
            self.check_entries(out, desired_entries)


# make a test suite to run all of the tests
loader = TestLoader()
suite = loader.loadTestsFromModule(sys.modules[__name__])


def test_all():
    TextTestRunner(verbosity=2).run(suite)

if __name__ == "__main__":
    test_all()
