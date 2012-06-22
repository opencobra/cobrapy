import unittest
import sys
from cobra.test import create_test_model
from cobra import Model, Reaction, Metabolite
from cobra.flux_analysis import *
from cobra.manipulation import initialize_growth_medium


class TestCobraFluxAnalysis(unittest.TestCase):
    """Test the simulation functions in cobra.flux_analysis

    TODO: Add in tests for: MOMA

    """
    def setUp(self):
        self.model = create_test_model()



    def test_single_deletion(self):
        cobra_model = self.model
        initialize_growth_medium(cobra_model, 'LB')

        #Expected growth rates for the salmonella model with deletions in LB medium
        the_loci =  ['STM4081', 'STM0247', 'STM3867', 'STM2952']
        the_genes = tpiA, metN, atpA, eno = map(cobra_model.genes.get_by_id, the_loci)
        id_to_name = dict([(x.id, x.name) for x in the_genes])
        growth_dict = {'fba':{tpiA.id:2.41, metN.id:2.44, atpA.id:1.87, eno.id:1.81}}


        for method, the_growth_rates in growth_dict.items():
            element_list = the_growth_rates.keys()
            rates, statuses, problems = single_deletion(cobra_model,
                                                        element_list=element_list,
                                                        method=method)

            for the_gene in element_list:
                self.assertEqual(statuses[the_gene], 'optimal')
                self.assertAlmostEqual(rates[the_gene], the_growth_rates[the_gene],
                                       places=2)


    def test_double_deletion(self):
        """
        """
        
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
        """

        """
        fva_results = {'5DGLCNtex': {'minimum': -1.9748300208638403e-05, 'maximum': 0.0}, 'ABTA': {'minimum': 0.0, 'maximum': 0.00014811225541408996}, '5DOAN': {'minimum': 0.0, 'maximum': 3.2227507421302166e-06}, 'A5PISO': {'minimum': 0.006920856282000001, 'maximum': 0.006922717378372606}, 'AACPS1': {'minimum': 0.0, 'maximum': 3.7028063376249126e-05}, 'AACPS2': {'minimum': 0.0, 'maximum': 3.7028063733878864e-05}, 'ACALDtex': {'minimum': -0.00011848980305159615, 'maximum': 0.0}, 'AACPS3': {'minimum': 0.0, 'maximum': 3.702806337623859e-05}, 'AACPS4': {'minimum': 0.0, 'maximum': 3.702806373387888e-05}, 'ABUTD': {'minimum': 0.0, 'maximum': 0.00014811225541406058}, 'AACPS5': {'minimum': 0.0, 'maximum': 2.8211857518389774e-05}, 'AACPS6': {'minimum': 0.0, 'maximum': 2.821185753295664e-05}, 'AACPS7': {'minimum': 0.0, 'maximum': 3.702806368868028e-05}, 'AACPS8': {'minimum': 0.0, 'maximum': 3.702806338788376e-05}, 'AACPS9': {'minimum': 0.0, 'maximum': 3.702806309933293e-05}, 'AACTOOR': {'minimum': 0.0, 'maximum': 1.5388286124597477e-05}, 'ABUTt2pp': {'minimum': 0.0, 'maximum': 0.0}, '3OAS140': {'minimum': 0.5041754136687804, 'maximum': 0.5042009621703677}, '3OAS141': {'minimum': 0.037484893950000084, 'maximum': 0.03750284695065363}, '3OAS160': {'minimum': 0.41767086529953557, 'maximum': 0.41769641380045963}, '3OAS161': {'minimum': 0.03748489395, 'maximum': 0.03750284695060761}, '3OAS180': {'minimum': 0.01069201669939239, 'maximum': 0.010717565200387778}, '3OAS181': {'minimum': 0.01606495455, 'maximum': 0.01608290755044158}, 'ABUTtex': {'minimum': 0.0, 'maximum': 0.0}, '3OAS60': {'minimum': 0.5439852127139995, 'maximum': 0.5439896193596934}, '3OAS80': {'minimum': 0.5439852127140001, 'maximum': 0.5439896193596934}, 'AAMYL': {'minimum': 0.0, 'maximum': 0.0}, '3PEPTabcpp': {'minimum': 0.0, 'maximum': 5.808323730923103e-06}, '3PEPTtex': {'minimum': -3.4245609402880297e-06, 'maximum': 0.0}, '3UMPtex': {'minimum': 0.0, 'maximum': 0.0}, '4HOXPACDtex': {'minimum': 0.0, 'maximum': 0.0}, 'ACACtex': {'minimum': 0.0, 'maximum': 0.0}, '4PCP': {'minimum': 0.0, 'maximum': 6.171343917391756e-06}, '4PCPpp': {'minimum': 0.0, 'maximum': 5.58914186256664e-06}, 'AAMYLpp': {'minimum': 0.0, 'maximum': 0.0}, '4PEPTabcpp': {'minimum': 0.0, 'maximum': 5.696625084349692e-06}, '4PEPTtex': {'minimum': -3.2198316806921494e-06, 'maximum': 0.0}, '5DGLCNR': {'minimum': -2.1942555793285538e-05, 'maximum': 0.0}, '5DGLCNt2rpp': {'minimum': -1.9748300208638403e-05, 'maximum': 0.0}, 'ACALD': {'minimum': 3.356574143593833, 'maximum': 7.4971939913624155}}


        cobra_model = self.model
        the_problem='return'
        initialize_growth_medium(cobra_model, 'LB')
        the_problem = cobra_model.optimize(the_problem=the_problem)
        fva_out = flux_variability_analysis(cobra_model,
                                            the_problem=the_problem,
                                            the_reactions=cobra_model.reactions[100:140])
        fva_out = dict([(k.id, v) for k, v in fva_out.items()])
        for the_reaction, the_range in fva_out.iteritems():
            for k, v in the_range.iteritems():
                self.assertAlmostEqual(fva_results[the_reaction][k], v, places=3)
        ## blocked_reactions = find_blocked_reactions(cobra_model)
        ## open_ex_blocked = find_blocked_reactions(cobra_model,
        ##                                          open_exchanges=True)

    def test_assess_medium_component_essentiality(self):
        """

        TODO: Add in a numerical value test
        """
        essentiality_results = {'EX_ser__L_e': 0.28511251509333996, 'EX_cobalt2_e': 0.0, 'EX_glu__L_e': 0.18551423955187463, 'EX_glyc_e': 0.02162967396132975, 'EX_h_e': 0.18551423955211313, 'EX_mobd_e': 0.0, 'EX_val__L_e': 0.18004717981556367, 'EX_so4_e': 0.1800471798156284, 'EX_co2_e': 0.18004717981574314, 'EX_k_e': 5.048709793414476e-29, 'EX_fe3_e': -3.4331226595218434e-27, 'EX_na1_e': 0.18004717981556367, 'EX_cl_e': 1.7495455763604752e-28, 'EX_leu__L_e': 0.1762785172191746, 'EX_arg__L_e': 0.13025755872698241, 'EX_nh4_e': 0.09432269135782297, 'EX_lys__L_e': 0.43718843672718055, 'EX_ala__L_e': 0.4371884367334397, 'EX_thr__L_e': 0.43718843673877533, 'EX_pi_e': 4.1028325973665373e-13, 'EX_mn2_e': 0.0, 'EX_phe__L_e': 0.380007972274807, 'EX_h2o_e': 0.38000797227380473, 'EX_mg2_e': 0.0, 'EX_his__L_e': 0.38000797227522415, 'EX_o2_e': 0.3428169207281707, 'EX_pro__L_e': 0.271070547843646, 'EX_asp__L_e': 0.38000797227507915, 'EX_gly_e': 0.3800079722747013, 'EX_cys__L_e': 0.3800079722760569, 'EX_cu2_e': 9.244463733058732e-31, 'EX_ca2_e': 0.0, 'EX_tyr__L_e': 0.38000797227331706, 'EX_zn2_e': 0.0, 'EX_met__L_e': 0.38000797227265026, 'EX_ile__L_e': 0.3800079722724871}
        cobra_model = self.model
        the_problem='return'
        from cobra.flux_analysis import assess_medium_component_essentiality
        essentiality_dict = assess_medium_component_essentiality(cobra_model, None, 'MgM')
        for k, v in essentiality_dict.iteritems():
            self.assertAlmostEqual(essentiality_results[k], v, places=3)




# make a test suite to run all of the tests
loader = unittest.TestLoader()
suite = loader.loadTestsFromModule(sys.modules[__name__])

def test_all():
    unittest.TextTestRunner(verbosity=2).run(suite)

if __name__ == "__main__":
    test_all()


