from unittest import TestCase, TestLoader, TextTestRunner, skipIf

from warnings import warn
import sys
from os.path import join
from json import load

from six import iteritems

try:
    import numpy
except:
    numpy = None

if __name__ == "__main__":
    sys.path.insert(0, "../..")
    from cobra.test import create_test_model, data_directory
    from cobra import Model, Reaction, Metabolite
    from cobra.manipulation import initialize_growth_medium
    from cobra.solvers import solver_dict, get_solver_name
    from cobra.manipulation import modify, delete
    from cobra.flux_analysis import *
    sys.path.pop(0)
else:
    from . import create_test_model, data_directory
    from .. import Model, Reaction, Metabolite
    from ..manipulation import initialize_growth_medium
    from ..solvers import solver_dict, get_solver_name
    from ..manipulation import modify, delete
    from ..flux_analysis import *


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

    def test_modify_reversible(self):
        model1 = create_test_model("textbook")
        model1.optimize()
        model2 = create_test_model("textbook")
        modify.convert_to_irreversible(model2)
        model2.optimize()
        self.assertAlmostEqual(model1.solution.f, model2.solution.f, places=3)
        modify.revert_to_reversible(model2)
        model2.optimize()
        self.assertAlmostEqual(model1.solution.f, model2.solution.f, places=3)

        # Ensure revert_to_reversible is robust to solutions generated both
        # before and after reversibility conversion, or not solved at all.
        model3 = create_test_model("textbook")
        model3.optimize()
        modify.convert_to_irreversible(model3)
        modify.revert_to_reversible(model3)
        self.assertAlmostEqual(model1.solution.f, model3.solution.f, places=3)

        model4 = create_test_model("textbook")
        modify.convert_to_irreversible(model4)
        modify.revert_to_reversible(model4)

    def test_escape_ids(self):
        model = create_test_model('textbook')
        model.reactions.PGI.gene_reaction_rule = "a.b or c"
        self.assertIn("a.b", model.genes)
        modify.escape_ID(model)
        self.assertNotIn("a.b", model.genes)

    def test_gene_knockout_computation(self):
        cobra_model = create_test_model()

        # helper functions for running tests
        delete_model_genes = delete.delete_model_genes
        find_gene_knockout_reactions = delete.find_gene_knockout_reactions

        def find_gene_knockout_reactions_fast(cobra_model, gene_list):
            compiled_rules = delete.get_compiled_gene_reaction_rules(
                cobra_model)
            return find_gene_knockout_reactions(
                cobra_model, gene_list,
                compiled_gene_reaction_rules=compiled_rules)

        def get_removed(m):
            return {x.id for x in m._trimmed_reactions}

        def test_computation(m, gene_ids, expected_reaction_ids):
            genes = [m.genes.get_by_id(i) for i in gene_ids]
            expected_reactions = {m.reactions.get_by_id(i)
                                  for i in expected_reaction_ids}
            removed1 = set(find_gene_knockout_reactions(m, genes))
            removed2 = set(find_gene_knockout_reactions_fast(m, genes))
            self.assertEqual(removed1, expected_reactions)
            self.assertEqual(removed2, expected_reactions)
            delete.delete_model_genes(m, gene_ids, cumulative_deletions=False)
            self.assertEqual(get_removed(m), expected_reaction_ids)
            delete.undelete_model_genes(m)

        gene_list = ['STM1067', 'STM0227']
        dependent_reactions = {'3HAD121', '3HAD160', '3HAD80', '3HAD140',
                               '3HAD180', '3HAD100', '3HAD181', '3HAD120',
                               '3HAD60', '3HAD141', '3HAD161', 'T2DECAI',
                               '3HAD40'}
        test_computation(cobra_model, gene_list, dependent_reactions)
        test_computation(cobra_model, ['STM4221'], {'PGI'})
        test_computation(cobra_model, ['STM1746.S'], {'4PEPTabcpp'})
        # test cumulative behavior
        delete_model_genes(cobra_model, gene_list[:1])
        delete_model_genes(cobra_model, gene_list[1:],
                           cumulative_deletions=True)
        delete_model_genes(cobra_model, ["STM4221"],
                           cumulative_deletions=True)
        dependent_reactions.add('PGI')
        self.assertEqual(get_removed(cobra_model), dependent_reactions)
        # non-cumulative following cumulative
        delete_model_genes(cobra_model, ["STM4221"],
                           cumulative_deletions=False)
        self.assertEqual(get_removed(cobra_model), {'PGI'})
        # make sure on reset that the bounds are correct
        reset_bound = cobra_model.reactions.get_by_id("T2DECAI").upper_bound
        self.assertEqual(reset_bound, 1000.)
        # test computation when gene name is a subset of another
        test_model = Model()
        test_reaction_1 = Reaction("test1")
        test_reaction_1.gene_reaction_rule = "eggs or (spam and eggspam)"
        test_model.add_reaction(test_reaction_1)
        test_computation(test_model, ["eggs"], set())
        test_computation(test_model, ["eggs", "spam"], {'test1'})
        # test computation with nested boolean expression
        test_reaction_1.gene_reaction_rule = \
            "g1 and g2 and (g3 or g4 or (g5 and g6))"
        test_computation(test_model, ["g3"], set())
        test_computation(test_model, ["g1"], {'test1'})
        test_computation(test_model, ["g5"], set())
        test_computation(test_model, ["g3", "g4", "g5"], {'test1'})
        # test computation when gene names are python expressions
        test_reaction_1.gene_reaction_rule = "g1 and (for or in)"
        test_computation(test_model, ["for", "in"], {'test1'})
        test_computation(test_model, ["for"], set())
        test_reaction_1.gene_reaction_rule = "g1 and g2 and g2.conjugate"
        test_computation(test_model, ["g2"], {"test1"})
        test_computation(test_model, ["g2.conjugate"], {"test1"})
        test_reaction_1.gene_reaction_rule = "g1 and (try:' or 'except:1)"
        test_computation(test_model, ["try:'"], set())
        test_computation(test_model, ["try:'", "'except:1"], {"test1"})

    def test_single_gene_deletion(self):
        cobra_model = create_test_model("textbook")
        # expected knockouts for textbook model
        growth_dict = {"fba": {"b0008": 0.87, "b0114": 0.80, "b0116": 0.78,
                               "b2276": 0.21, "b1779": 0.00},
                       "moma": {"b0008": 0.87, "b0114": 0.71, "b0116": 0.56,
                                "b2276": 0.11, "b1779": 0.00},
                       }

        # MOMA requires cplex or gurobi
        try:
            get_solver_name(qp=True)
        except:
            growth_dict.pop('moma')
        for method, expected in growth_dict.items():
            rates, statuses = single_gene_deletion(cobra_model,
                                                   gene_list=expected.keys(),
                                                   method=method)
            for gene, expected_value in iteritems(expected):
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
        solution = double_gene_deletion(cobra_model, gene_list1=genes)
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
            self.skip("no MILP solver found")
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
            self.skip("no MILP solver found")
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

        result = gapfilling.growMatch(m, U)[0]
        self.assertEqual(len(result), 1)
        self.assertEqual(result[0].id, "a2b")

        # 2 rounds with exchange reactions
        result = gapfilling.growMatch(m, None, ex_rxns=True, iterations=2)
        self.assertEqual(len(result), 2)
        self.assertEqual(len(result[0]), 1)
        self.assertEqual(len(result[1]), 1)
        self.assertEqual({i[0].id for i in result},
                         {"SMILEY_EX_b", "SMILEY_EX_c"})

    def test_remove_genes(self):
        m = Model("test")
        m.add_reactions([Reaction("r" + str(i + 1)) for i in range(8)])
        self.assertEqual(len(m.reactions), 8)
        rxns = m.reactions
        rxns.r1.gene_reaction_rule = "(a and b) or (c and a)"
        rxns.r2.gene_reaction_rule = "(a and b and d and e)"
        rxns.r3.gene_reaction_rule = "(a and b) or (b and c)"
        rxns.r4.gene_reaction_rule = "(f and b) or (b and c)"
        rxns.r5.gene_reaction_rule = "x"
        rxns.r6.gene_reaction_rule = "y"
        rxns.r7.gene_reaction_rule = "x or     z"
        rxns.r8.gene_reaction_rule = ""
        self.assertIn("a", m.genes)
        self.assertIn("x", m.genes)
        delete.remove_genes(m, ["a"], remove_reactions=False)
        self.assertNotIn("a", m.genes)
        self.assertIn("x", m.genes)
        self.assertEqual(rxns.r1.gene_reaction_rule, "")
        self.assertEqual(rxns.r2.gene_reaction_rule, "")
        self.assertEqual(rxns.r3.gene_reaction_rule, "b and c")
        self.assertEqual(rxns.r4.gene_reaction_rule, "(f and b) or (b and c)")
        self.assertEqual(rxns.r5.gene_reaction_rule, "x")
        self.assertEqual(rxns.r6.gene_reaction_rule, "y")
        self.assertEqual(rxns.r7.genes, {m.genes.x, m.genes.z})
        self.assertEqual(rxns.r8.gene_reaction_rule, "")
        delete.remove_genes(m, ["x"], remove_reactions=True)
        self.assertEqual(len(m.reactions), 7)
        self.assertNotIn("r5", m.reactions)
        self.assertNotIn("x", m.genes)
        self.assertEqual(rxns.r1.gene_reaction_rule, "")
        self.assertEqual(rxns.r2.gene_reaction_rule, "")
        self.assertEqual(rxns.r3.gene_reaction_rule, "b and c")
        self.assertEqual(rxns.r4.gene_reaction_rule, "(f and b) or (b and c)")
        self.assertEqual(rxns.r6.gene_reaction_rule, "y")
        self.assertEqual(rxns.r7.gene_reaction_rule, "z")
        self.assertEqual(rxns.r7.genes, {m.genes.z})
        self.assertEqual(rxns.r8.gene_reaction_rule, "")

    @skipIf(numpy is None, "double deletions require numpy")
    def test_phenotype_phase_plane(self):
        model = create_test_model("textbook")
        data = calculate_phenotype_phase_plane(
            model, "EX_glc__D_e", "EX_o2_e",
            reaction1_npoints=20, reaction2_npoints=20)
        self.assertEqual(data.growth_rates.shape, (20, 20))
        self.assertAlmostEqual(data.growth_rates.max(), 1.20898, places=4)
        self.assertAlmostEqual(abs(data.growth_rates[0, :]).max(), 0, places=4)

# make a test suite to run all of the tests
loader = TestLoader()
suite = loader.loadTestsFromModule(sys.modules[__name__])


def test_all():
    TextTestRunner(verbosity=2).run(suite)

if __name__ == "__main__":
    test_all()
