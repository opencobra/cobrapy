from unittest import TestCase, TestLoader, TextTestRunner, skipIf

from warnings import warn
import sys
from os import name

if __name__ == "__main__":
    sys.path.insert(0, "../..")
    from cobra.test import create_test_model
    from cobra import Model, Reaction, Metabolite
    from cobra.manipulation import initialize_growth_medium
    from cobra.solvers import solver_dict
    from cobra.manipulation import modify
    from cobra.flux_analysis.parsimonious import optimize_minimal_flux
    from cobra.flux_analysis.single_deletion import single_deletion
    from cobra.flux_analysis.variability import flux_variability_analysis
    from cobra.flux_analysis.variability import flux_variability_analysis_fast
    sys.path.pop(0)
else:
    from . import create_test_model
    from .. import Model, Reaction, Metabolite
    from ..manipulation import initialize_growth_medium
    from ..solvers import solver_dict
    from ..manipulation import modify
    from ..flux_analysis.parsimonious import optimize_minimal_flux
    from ..flux_analysis.single_deletion import single_deletion
    from ..flux_analysis.double_deletion import double_deletion
    from ..flux_analysis.variability import flux_variability_analysis
    from ..flux_analysis.variability import flux_variability_analysis_fast


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
            self.assertGreater(min(abs_x), -1e-7)

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


    def test_single_deletion(self):
        cobra_model = self.model
        initialize_growth_medium(cobra_model, 'LB')

        #Expected growth rates for the salmonella model with deletions in LB medium
        the_loci =  ['STM4081', 'STM0247', 'STM3867', 'STM2952']
        the_genes = tpiA, metN, atpA, eno = map(cobra_model.genes.get_by_id, the_loci)
        id_to_name = dict([(x.id, x.name) for x in the_genes])
        growth_dict = {'fba':{tpiA.id:2.41, metN.id:2.44, atpA.id:1.87, eno.id:1.81},
                       'moma':{ tpiA.id:1.62, metN.id:2.4, atpA.id:1.40, eno.id:0.33}}

        #MOMA requires cplex or gurobi
        if 'cplex' not in solver_dict and 'gurobi' not in solver_dict:
            growth_dict.pop('moma')
        for method, the_growth_rates in growth_dict.items():
            element_list = the_growth_rates.keys()
            results = single_deletion(cobra_model, element_list=element_list,
                                      method=method)
            rates = results[0]
            statuses = results[1]

            for the_gene in element_list:
                self.assertEqual(statuses[the_gene], 'optimal')
                self.assertAlmostEqual(rates[the_gene], the_growth_rates[the_gene],
                                       places=2)

    @skipIf(name == "java", "cobra.test.flux_analysis.test_double_deletion doesn't yet work with java")
    def test_double_deletion(self):
        cobra_model = self.model
        #turn into a double deletion unit test
        the_problem='return'
        initialize_growth_medium(cobra_model, 'LB')
        #Expected growth rates for the salmonella model with deletions in LB medium
        the_loci =  ['STM4081', 'STM0247', 'STM3867', 'STM2952']
        the_genes = tpiA, metN, atpA, eno = map(cobra_model.genes.get_by_id, the_loci)
        growth_dict = {}
        growth_list = [[2.41, 2.389, 1.775, 1.81],
                       [2.389, 2.437, 1.86, 1.79],
                       [1.775, 1.86, 1.87, 1.3269],
                       [1.81, 1.79, 1.3269, 1.81]]
        for the_gene, the_rates in zip(the_genes, growth_list):
            growth_dict[the_gene] = dict(zip(the_genes, the_rates))


        the_solution = double_deletion(cobra_model, element_list_1=the_genes,
                                       element_list_2=the_genes,
                                       the_problem=the_problem)
        #Potential problem if the data object doesn't have a tolist function
        s_data = the_solution['data'].tolist()
        s_x = the_solution['x']
        s_y = the_solution['y']
        for gene_x, rates_x in zip(s_x, s_data):
            for gene_y, the_rate in zip(s_y, rates_x):
                self.assertAlmostEqual(growth_dict[gene_x][gene_y], the_rate,
                                       places=2)

    def test_flux_variability(self):
        fva_results = {
            '5DGLCNtex': {'minimum': -1.9748300208638403e-05, 'maximum': 0.0},
            'ABTA': {'minimum': 0.0, 'maximum': 0.00014811225541408996},
            '5DOAN': {'minimum': 0.0, 'maximum': 3.2227507421302166e-06},
            'A5PISO': {'minimum': 0.006920856282000001, 'maximum': 0.006922717378372606},
            'AACPS1': {'minimum': 0.0, 'maximum': 3.7028063376249126e-05},
            'AACPS2': {'minimum': 0.0, 'maximum': 3.7028063733878864e-05},
            'ACALDtex': {'minimum': -0.00011848980305159615, 'maximum': 0.0},
            'AACPS3': {'minimum': 0.0, 'maximum': 3.702806337623859e-05},
            'AACPS4': {'minimum': 0.0, 'maximum': 3.702806373387888e-05},
            'ABUTD': {'minimum': 0.0, 'maximum': 0.00014811225541406058},
            'AACPS5': {'minimum': 0.0, 'maximum': 2.8211857518389774e-05},
            'AACPS6': {'minimum': 0.0, 'maximum': 2.821185753295664e-05},
            'AACPS7': {'minimum': 0.0, 'maximum': 3.702806368868028e-05},
            'AACPS8': {'minimum': 0.0, 'maximum': 3.702806338788376e-05},
            'AACPS9': {'minimum': 0.0, 'maximum': 3.702806309933293e-05},
            'AACTOOR': {'minimum': 0.0, 'maximum': 1.5388286124597477e-05},
            'ABUTt2pp': {'minimum': 0.0, 'maximum': 0.0},
            '3OAS140': {'minimum': 0.5041754136687804, 'maximum': 0.5042009621703677},
            '3OAS141': {'minimum': 0.037484893950000084, 'maximum': 0.03750284695065363},
            '3OAS160': {'minimum': 0.41767086529953557, 'maximum': 0.41769641380045963},
            '3OAS161': {'minimum': 0.03748489395, 'maximum': 0.03750284695060761},
            '3OAS180': {'minimum': 0.01069201669939239, 'maximum': 0.010717565200387778},
            '3OAS181': {'minimum': 0.01606495455, 'maximum': 0.01608290755044158},
            'ABUTtex': {'minimum': 0.0, 'maximum': 0.0},
            '3OAS60': {'minimum': 0.5439852127139995, 'maximum': 0.5439896193596934},
            '3OAS80': {'minimum': 0.5439852127140001, 'maximum': 0.5439896193596934},
            'AAMYL': {'minimum': 0.0, 'maximum': 0.0},
            '3PEPTabcpp': {'minimum': 0.0, 'maximum': 5.808323730923103e-06},
            '3PEPTtex': {'minimum': -3.4245609402880297e-06, 'maximum': 0.0},
            '3UMPtex': {'minimum': 0.0, 'maximum': 0.0},
            '4HOXPACDtex': {'minimum': 0.0, 'maximum': 0.0},
            'ACACtex': {'minimum': 0.0, 'maximum': 0.0},
            '4PCP': {'minimum': 0.0, 'maximum': 6.171343917391756e-06},
            '4PCPpp': {'minimum': 0.0, 'maximum': 5.58914186256664e-06},
            'AAMYLpp': {'minimum': 0.0, 'maximum': 0.0},
            '4PEPTabcpp': {'minimum': 0.0, 'maximum': 5.696625084349692e-06},
            '4PEPTtex': {'minimum': -3.2198316806921494e-06, 'maximum': 0.0},
            '5DGLCNR': {'minimum': -2.1942555793285538e-05, 'maximum': 0.0},
            '5DGLCNt2rpp': {'minimum': -1.9748300208638403e-05, 'maximum': 0.0},
            'ACALD': {'minimum': 3.356574143593833, 'maximum': 7.4971939913624155}}

        cobra_model = self.model
        the_problem='return'
        initialize_growth_medium(cobra_model, 'LB')
        the_problem = cobra_model.optimize(the_problem=the_problem)
        fva_out = flux_variability_analysis(cobra_model,
                                            the_problem=the_problem,
                                            the_reactions=cobra_model.reactions[100:140])
        for the_reaction, the_range in fva_out.iteritems():
            for k, v in the_range.iteritems():
                self.assertAlmostEqual(fva_results[the_reaction][k], v, places=3)
        # test fastfva - only works to 2 decimal places
        fva_out = flux_variability_analysis_fast(cobra_model,
                                                 cobra_model.reactions[100:140],
                                                 fraction_of_optimum=1.)
        for the_reaction, the_range in fva_out.iteritems():
            for k, v in the_range.iteritems():
                self.assertAlmostEqual(fva_results[the_reaction][k], v, places=2)

# make a test suite to run all of the tests
loader = TestLoader()
suite = loader.loadTestsFromModule(sys.modules[__name__])

def test_all():
    TextTestRunner(verbosity=2).run(suite)

if __name__ == "__main__":
    test_all()
