from unittest import TestCase, TestLoader, TextTestRunner, skipIf

from warnings import warn
import sys
from os import name

try:
    import numpy
except:
    numpy = None

if __name__ == "__main__":
    sys.path.insert(0, "../..")
    from cobra.test import create_test_model
    from cobra import Model, Reaction, Metabolite
    from cobra.manipulation import initialize_growth_medium
    from cobra.solvers import solver_dict, get_solver_name
    from cobra.manipulation import modify, delete
    from cobra.flux_analysis.parsimonious import optimize_minimal_flux
    from cobra.flux_analysis.variability import flux_variability_analysis
    from cobra.flux_analysis.variability import find_blocked_reactions
    from cobra.flux_analysis.single_deletion import single_deletion
    from cobra.flux_analysis.loopless import construct_loopless_model
    from cobra.external.six import iteritems
    if numpy:
        from cobra.flux_analysis.double_deletion import double_deletion
    sys.path.pop(0)
else:
    from . import create_test_model
    from .. import Model, Reaction, Metabolite
    from ..manipulation import initialize_growth_medium
    from ..solvers import solver_dict, get_solver_name
    from ..manipulation import modify, delete
    from ..flux_analysis.parsimonious import optimize_minimal_flux
    from ..flux_analysis.variability import flux_variability_analysis
    from ..flux_analysis.variability import find_blocked_reactions
    from ..flux_analysis.single_deletion import single_deletion
    from ..flux_analysis.loopless import construct_loopless_model
    from ..external.six import iteritems
    if numpy:
        from ..flux_analysis.double_deletion import double_deletion


class TestCobraFluxAnalysis(TestCase):
    """Test the simulation functions in cobra.flux_analysis"""

    def setUp(self):
        self.model = create_test_model()

    def test_pFBA(self):
        model = self.model
        for solver in solver_dict:
            optimize_minimal_flux(model, solver=solver)
            abs_x = [abs(i) for i in model.solution.x]
            self.assertAlmostEqual(model.solution.f, 0.3800, places=3)
            self.assertAlmostEqual(sum(abs_x), 343.021, places=3)

    def test_modify_reversible(self):
        model1 = self.model
        model1.optimize()
        model2 = create_test_model()
        modify.convert_to_irreversible(model2)
        model2.optimize()
        self.assertAlmostEqual(model1.solution.f, model2.solution.f, places=3)
        modify.revert_to_reversible(model2)
        model2.optimize()
        self.assertAlmostEqual(model1.solution.f, model2.solution.f, places=3)

        # Ensure revert_to_reversible is robust to solutions generated both
        # before and after reversibility conversion, or not solved at all.
        model3 = create_test_model()
        model3.optimize()
        modify.convert_to_irreversible(model3)
        modify.revert_to_reversible(model3)
        self.assertAlmostEqual(model1.solution.f, model3.solution.f, places=3)

        model4 = create_test_model()
        modify.convert_to_irreversible(model4)
        modify.revert_to_reversible(model4)

    def test_gene_knockout_computation(self):
        cobra_model = self.model

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

    def test_single_deletion(self):
        cobra_model = self.model
        initialize_growth_medium(cobra_model, 'LB')

        # Expected growth rates for the salmonella model with deletions in LB
        the_loci = ['STM4081', 'STM0247', 'STM3867', 'STM2952']
        the_genes = tpiA, metN, atpA, eno = map(cobra_model.genes.get_by_id,
                                                the_loci)
        id_to_name = dict([(x.id, x.name) for x in the_genes])
        growth_dict = {'fba': {tpiA.id: 2.41, metN.id: 2.44,
                               atpA.id: 1.87, eno.id: 1.81},
                       'moma': {tpiA.id: 1.62, metN.id: 2.4,
                                atpA.id: 1.40, eno.id: 0.33}}

        # MOMA requires cplex or gurobi
        try:
            get_solver_name(qp=True)
        except:
            growth_dict.pop('moma')
        for method, the_growth_rates in growth_dict.items():
            element_list = the_growth_rates.keys()
            results = single_deletion(cobra_model, element_list=element_list,
                                      method=method)
            rates = results[0]
            statuses = results[1]

            for the_gene in element_list:
                self.assertEqual(statuses[the_gene], 'optimal')
                self.assertAlmostEqual(rates[the_gene],
                                       the_growth_rates[the_gene],
                                       places=2)

    def test_single_reaction_deletion(self):
        cobra_model = self.model
        results, status = single_deletion(
            cobra_model, element_list=cobra_model.reactions[100:110:2],
            element_type="reaction")
        self.assertEqual(len(results), 5)
        self.assertEqual(len(status), 5)
        for status_value in status.values():
            self.assertEqual(status_value, "optimal")
        expected_results = {'3OAS140': 0, '3OAS160': 0.38001,
                            '3OAS180': 0.38001, '3OAS60': 0,
                            '3PEPTabcpp': 0.38001}
        for reaction, value in results.items():
            self.assertAlmostEqual(value, expected_results[reaction], 5)

    @skipIf(numpy is None, "double deletions require numpy")
    def test_double_deletion(self):
        cobra_model = self.model
        # turn into a double deletion unit test
        initialize_growth_medium(cobra_model, 'LB')
        # Expected growth rates for the salmonella model with deletions in LB
        the_loci = ['STM4081', 'STM0247', 'STM3867', 'STM2952']
        the_genes = list(map(cobra_model.genes.get_by_id, the_loci))
        growth_dict = {}
        growth_list = [[2.41, 2.389, 1.775, 1.81],
                       [2.389, 2.437, 1.86, 1.79],
                       [1.775, 1.86, 1.87, 1.3269],
                       [1.81, 1.79, 1.3269, 1.81]]
        for the_gene, the_rates in zip(the_genes, growth_list):
            growth_dict[the_gene] = dict(zip(the_genes, the_rates))

        the_solution = double_deletion(cobra_model, element_list_1=the_genes,
                                       element_list_2=the_genes)
        # Potential problem if the data object doesn't have a tolist function
        s_data = the_solution['data'].tolist()
        s_x = the_solution['x']
        s_y = the_solution['y']
        for gene_x, rates_x in zip(s_x, s_data):
            for gene_y, the_rate in zip(s_y, rates_x):
                self.assertAlmostEqual(
                    growth_dict[gene_x][gene_y], the_rate, places=2,
                    msg="%.3f != %.3f for (%s, %s)" %
                    (growth_dict[gene_x][gene_y], the_rate,
                     gene_x.id, gene_y.id))

    def test_flux_variability(self):
        fva_results = {
            '5DGLCNtex': {'minimum': 0.0, 'maximum': 0.0},
            'ABTA': {'minimum': 0.0, 'maximum': 0.0},
            '5DOAN': {'minimum': 0.0, 'maximum': 0.0},
            'A5PISO': {'minimum': 0.00692, 'maximum': 0.00692},
            'AACPS1': {'minimum': 0.0, 'maximum': 0.0},
            'AACPS2': {'minimum': 0.0, 'maximum': 0.0},
            'ACALDtex': {'minimum': 0.0, 'maximum': 0.0},
            'AACPS3': {'minimum': 0.0, 'maximum': 0.0},
            'AACPS4': {'minimum': 0.0, 'maximum': 0.0},
            'ABUTD': {'minimum': 0.0, 'maximum': 0.0},
            'AACPS5': {'minimum': 0.0, 'maximum': 0.0},
            'AACPS6': {'minimum': 0.0, 'maximum': 0.0},
            'AACPS7': {'minimum': 0.0, 'maximum': 0.0},
            'AACPS8': {'minimum': 0.0, 'maximum': 0.0},
            'AACPS9': {'minimum': 0.0, 'maximum': 0.0},
            'AACTOOR': {'minimum': 0.0, 'maximum': 0.0},
            'ABUTt2pp': {'minimum': 0.0, 'maximum': 0.0},
            '3OAS140': {'minimum': 0.50419, 'maximum': 0.50419},
            '3OAS141': {'minimum': 0.03748, 'maximum': 0.03748},
            '3OAS160': {'minimum': 0.41769, 'maximum': 0.41769},
            '3OAS161': {'minimum': 0.03748, 'maximum': 0.03748},
            '3OAS180': {'minimum': 0.01071, 'maximum': 0.01071},
            '3OAS181': {'minimum': 0.01606, 'maximum': 0.01606},
            'ABUTtex': {'minimum': 0.0, 'maximum': 0.0},
            '3OAS60': {'minimum': 0.54399, 'maximum': 0.54399},
            '3OAS80': {'minimum': 0.54399, 'maximum': 0.54399},
            'AAMYL': {'minimum': 0.0, 'maximum': 0.0},
            '3PEPTabcpp': {'minimum': 0.0, 'maximum': 0.0},
            '3PEPTtex': {'minimum': 0.0, 'maximum': 0.0},
            '3UMPtex': {'minimum': 0.0, 'maximum': 0.0},
            '4HOXPACDtex': {'minimum': 0.0, 'maximum': 0.0},
            'ACACtex': {'minimum': 0.0, 'maximum': 0.0},
            '4PCP': {'minimum': 0.0, 'maximum': 0.0},
            '4PCPpp': {'minimum': 0.0, 'maximum': 0.0},
            'AAMYLpp': {'minimum': 0.0, 'maximum': 0.0},
            '4PEPTabcpp': {'minimum': 0.0, 'maximum': 0.0},
            '4PEPTtex': {'minimum': 0.0, 'maximum': 0.0},
            '5DGLCNR': {'minimum': 0.0, 'maximum': 0.0},
            '5DGLCNt2rpp': {'minimum': 0.0, 'maximum': 0.0},
            'ACALD': {'minimum': 3.35702, 'maximum': 7.49572}}

        infeasible_model = create_test_model()
        infeasible_model.reactions.get_by_id("EX_glyc_e").lower_bound = 0
        for solver in solver_dict:
            cobra_model = create_test_model()
            initialize_growth_medium(cobra_model, 'LB')
            fva_out = flux_variability_analysis(
                cobra_model, solver=solver,
                reaction_list=cobra_model.reactions[100:140:2])
            for the_reaction, the_range in iteritems(fva_out):
                for k, v in iteritems(the_range):
                    self.assertAlmostEqual(fva_results[the_reaction][k], v,
                                           places=5)
            # ensure that an infeasible model does not run FVA
            self.assertRaises(ValueError, flux_variability_analysis,
                              infeasible_model, solver=solver)

    def test_find_blocked_reactions(self):
        m = self.model
        result = find_blocked_reactions(m, reaction_list=m.reactions[:10])
        self.assertEqual(len(result), 1)
        self.assertEqual(result[0], '12PPDStex')

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

# make a test suite to run all of the tests
loader = TestLoader()
suite = loader.loadTestsFromModule(sys.modules[__name__])


def test_all():
    TextTestRunner(verbosity=2).run(suite)

if __name__ == "__main__":
    test_all()
