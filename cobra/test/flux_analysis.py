import unittest
import warnings
import os
try:
    from cPickle import load
except:
    from pickle import load
import sys
raise Exception("flux_analysis unit test has not been implemented")

# from . import data_directory, ecoli_sbml, ecoli_pickle, create_test_model
from . import create_test_model
from .. import Model, Reaction, Metabolite
from .. import solvers


class TestCobraSolver(unittest.TestCase):
    def setUp(self):
        self.model = create_test_model()
        self.old_solution = 0.982371812727
        self.infeasible_problem = Model()
        metabolite_1 = Metabolite("met1")
        metabolite_2 = Metabolite("met2")
        reaction_1 = Reaction("rxn1")
        reaction_2 = Reaction("rxn2")
        reaction_1.add_metabolites({metabolite_1: 1})
        reaction_2.add_metabolites({metabolite_2: 1})
        reaction_1.lower_bound = 1
        reaction_2.upper_bound = 2
        self.infeasible_problem.add_reactions([reaction_1, reaction_2])
        self.infeasible_problem.update()


def add_test(TestCobraSolver, solver_name, solver):
    def test_attributes(self):
        self.assertTrue(hasattr(solver, "create_problem"))
        self.assertTrue(hasattr(solver, "solve_problem"))
        # self.assertTrue(hasattr(solver, "update_problem"))
    def test_setup(self):
        solver.create_problem(self.model)
    def test_solve_feasible(self):
        solution = solver.solve(self.model)
        self.assertEqual(solution.status, "optimal")
        self.assertAlmostEqual(self.old_solution, \
            solution.objective_value, places=4)
    def test_solve_infeasible(self):
        solution = solver.solve(self.infeasible_problem)
        self.assertEqual(solution.status, "infeasible")
    def test_independent_creation(self):
        feasible_lp = solver.create_problem(self.model)
        infeasible_lp = solver.create_problem(self.infeasible_problem)
        feasible_solution = solve_problem(lp)
        infeasible_solution = solve_problem(lp)
        self.assertEqual(feasible_solution.status, "optimal")
        self.assertAlmostEqual(self.old_solution, \
            feasible_solution.objective_value, places=4)
        self.assertEqual(infeasible_solution.status, "infeasible")
    setattr(TestCobraSolver, "test_%s_create" % solver_name, \
        test_setup)
    setattr(TestCobraSolver, "test_%s_attributes" % solver_name, \
        test_attributes)
    setattr(TestCobraSolver, "test_%s_feasible_solve" % solver_name, \
        test_solve_feasible)
    setattr(TestCobraSolver, "test_%s_infeasible_solve" % solver_name, \
        test_solve_infeasible)
    setattr(TestCobraSolver, "test_%s_independent_creation" % solver_name, \
        test_solve_infeasible)

for solver_name, solver in solvers.solver_list.iteritems():
    add_test(TestCobraSolver, solver_name, solver)
# make a test suite to run all of the tests
loader = unittest.TestLoader()
suite = loader.loadTestsFromModule(sys.modules[__name__])

def test_all():
    unittest.TextTestRunner(verbosity=2).run(suite)

if __name__ == "__main__":
    test_all()


######################
#This goes in the variablity analysis test
if __name__ == '__main__':
    from cPickle import load
    from time import time
    from cobra.test import salmonella_pickle
    with open(salmonella_pickle) as in_file:
        cobra_model = load(in_file)


    solver_dict = {'glpk': True,
                   'gurobi': True,
                   'cplex': True}
    try:
        import glpk
    except:
        solver_dict['glpk'] = False
    try:
        from gurobipy import Model
    except:
        solver_dict['gurobi'] = False
    try:
        from cplex import Cplex
    except:
        solver_dict['cplex'] = False
    #The next two lines are used when profiling / debugging
    #from cProfile import run
    #for solver in solver_dict.items()[:0]:
    for solver, solver_status in solver_dict.items():
        if solver_status:
            print '\ntesting %s:'%solver
            print '\tFlux Variability: '
            the_problem = cobra_model.optimize(the_problem='return', solver=solver)
            start_time = time()
            fva_out = flux_variability_analysis(cobra_model,
                                                the_problem=the_problem,
                                                the_reactions=cobra_model.reactions,#[100:140],
                                                solver=solver,number_of_processes=1)
            print '\thot start: %f'%(time() - start_time)
            start_time = time()
            number_of_processes = 5
            pfva_out = flux_variability_analysis(cobra_model,
                                                the_problem=the_problem,
                                                the_reactions=cobra_model.reactions,#[100:140],
                                                solver=solver,number_of_processes=number_of_processes)
            print '\t%i processes hot start: %f'%(number_of_processes, time() - start_time)
            print '\n\tFind Blocked Reactions:'
            start_time = time()
            blocked_reactions = find_blocked_reactions(cobra_model,
                                                       solver=solver)
            print '\t\t:Basic Model: %f'%(time() - start_time)
            start_time = time()
            open_ex_blocked = find_blocked_reactions(cobra_model,
                                                     open_exchanges=True)
            print '\t\t:Opened Exchanges: %f'%(time() - start_time)



######################
#This goes in the single deletion test
if __name__ == '__main__':
    from cPickle import load
    from time import time
    from math import floor
    from cobra.test import salmonella_pickle
    method='moma'
    the_problem='return'
    element_type='gene'
    error_reporting=None
    from cobra.manipulation import initialize_growth_medium
    with open(salmonella_pickle) as in_file:
        cobra_model = load(in_file)


    initialize_growth_medium(cobra_model, 'LB')
    #Expected growth rates for the salmonella model with deletions in LB medium
    the_loci =  ['STM4081', 'STM0247', 'STM3867', 'STM2952']
    the_genes = tpiA, metN, atpA, eno = map(cobra_model.genes.get_by_id, the_loci)
    id_to_name = dict([(x.id, x.name) for x in the_genes])
    growth_dict = {'moma': {tpiA.id:1.61, metN.id:2.39, atpA.id:1.40, eno.id:0.33},
                   'fba':{tpiA.id:2.41, metN.id:2.43, atpA.id:1.87, eno.id:1.81}}

    solver_list = ['glpk',
                   'gurobi',
                   'cplex']
    try:
        import glpk
    except:
        solver_list.remove('glpk')
    try:
        from gurobipy import Model
    except:
        solver_list.remove('gurobi')
    try:
        from cplex import Cplex
    except:
        solver_list.remove('cplex')

    for solver in solver_list:
        print 'testing solver: ' + solver
        for method, the_growth_rates in growth_dict.items():
            if method == 'moma':
                print "Can't do MOMA now.  Try back later."
                continue
            print '\twith method: ' + method 
            element_list = the_growth_rates.keys()
            start_time = time()
            rates, statuses, problems = single_deletion(cobra_model,
                                                        element_list=element_list,
                                                        method=method,
                                                        the_problem=the_problem,
                                                        element_type=element_type,
                                                        solver=solver,
                                                        error_reporting=error_reporting)
            for the_gene, v in statuses.items():
                if v != 'optimal':
                    print '\t\tdeletion %s was not optimal'%the_gene
            for the_gene, v in rates.items():
                v = floor(100*v)/100
                if v != the_growth_rates[the_gene]:
                    print '\t\tFAILED: %s simulation (%1.3f) != expectation (%1.2f)'%(id_to_name[the_gene],
                                                                                     v,
                                                                                     the_growth_rates[the_gene])
                else:
                    print '\t\tPASSED: %s simulation (%1.3f) ~= expectation (%1.2f)'%(id_to_name[the_gene],
                                                                                     v,
                                                                                     the_growth_rates[the_gene])


            print '\t\tsingle deletion time: %f seconds'%(time() - start_time)




####################
#Merge this into the single and double deletion tests

if __name__ == '__main__':
    from cPickle import load
    from time import time
    from math import floor
    from numpy import round
    from cobra.flux_analysis.single_deletion import single_deletion
    from cobra.flux_analysis.double_deletion import double_deletion
    from cobra.test import salmonella_pickle
    objective_sense='minimize'
    tolerance_optimality=1e-8
    tolerance_feasibility=1e-8
    minimize_norm=False
    relax_b=False
    print_time = True
    the_problem='return'
    combined_model=None
    norm_type = 'euclidean'
    try:
        from gurobipy import Model
        solver, lp_method = 'gurobi', 0
        print "Using solver: %s with lp_method: %i"%(solver, lp_method)
    except:
        print "Gurobi isn't available trying cplex"
        try:
            from cplex import Cplex
            solver, lp_method = 'cplex', 1
            print "Using solver: %s with lp_method: %i"%(solver, lp_method)
        except:
            print "Couldn't load gurobi or cplex so can't run moma trying lmoma"
            try:
                from glpk import LPX
                solver, lp_method = 'glpk', 1
                norm_type = 'linear'
                print "Using solver: %s with lp_method: %i"%(solver, lp_method)
            except:
                print "Couldn't load gurobi, cplex, or glpk"
                
    tpiA_result = 1.61
    tpiA_metN_result = 1.60
    with open(salmonella_pickle) as in_file:
        cobra_model = load(in_file)
    gene_list = ['tpiA', 'metN']
    loci_list =  ['STM4081', 'STM0247']
    initialize_growth_medium(cobra_model, 'LB')
    wt_model = cobra_model
    #gene_list = map(wt_model.genes.get_by_id, loci_list)
    wt_model.id = 'Wild-type'
    start_time = time()
    mutant_model = wt_model.copy()
    print 'copy time: %1.2f'%(time() - start_time)
    mutant_model.id = 'mutant %s'%gene_list[0]
    delete_model_genes(mutant_model, gene_list[:1])

    if norm_type == 'linear':
        start_time = time()
        the_solution = moma(wt_model, mutant_model, solver=solver,
                            lp_method=lp_method, norm_type=norm_type,
                            print_time=print_time)
        print time() - start_time
        the_problem = the_solution['the_problem']
        combined_model = the_solution['combined_model']
        print 'tpiA: %1.2f'%the_solution['objective_value']
        start_time = time()
        the_solution = moma(wt_model, mutant_model, solver=solver,
                            lp_method=lp_method, the_problem=the_problem,
                            norm_type=norm_type, print_time=print_time,
                            combined_model=combined_model)
        print 'tpiA: %1.2f'%the_solution['objective_value']
        print time() - start_time
        
    else:
        start_time = time()
        the_solution = moma(wt_model, mutant_model, solver=solver,
                            lp_method=lp_method, norm_type=norm_type,
                            print_time=print_time)
        if the_solution['status'] in ['optimal']:
            the_problem = the_solution['the_problem']
            tmp_result = floor(100*the_solution['objective_value'])/100
            if tmp_result == tpiA_result:
                print 'Passed MOMA with tpiA deletion in %1.4f seconds'%(time() - start_time) +\
                      '\n%1.2f == %1.2f simulated'%(tpiA_result,
                                                  tmp_result)
                
                
            else:
                print 'FAILED: tpiA deletion expected '+\
                      '%1.2f != %1.2f simulated'%(tpiA_result,
                                                  tmp_result)

                start_time = time()
                the_solution = moma(wt_model, mutant_model,
                                    the_problem=the_problem,
                                    solver=solver, lp_method=lp_method,
                                    norm_type=norm_type, print_time=print_time)
                if the_solution['status'] in ['optimal']:
                    the_problem = the_solution['the_problem']
                    tmp_result = floor(100*the_solution['objective_value'])/100
                    if tmp_result == tpiA_result:
                        print 'Passed MOMA reusing Model with tpiA deletion in %1.4f seconds'%(time() - start_time)
                start_time = time()
                the_solution = moma(wt_model, mutant_model,
                                    the_problem=the_problem,
                                    combined_model=the_solution['combined_model'], lp_method=lp_method, solver=solver,
                                    norm_type=norm_type, print_time=print_time)
                if the_solution['status'] in ['optimal']:
                    the_problem = the_solution['the_problem']
                    tmp_result = floor(100*the_solution['objective_value'])/100
                    if tmp_result == tpiA_result:
                        print 'Passed MOMA reusing Model and model with tpiA deletion in %1.4f seconds'%(time() - start_time)
            start_time - time()
            single_solution = single_deletion(wt_model, loci_list, method='moma', solver=solver)
            the_status = single_solution[1]
            if the_status['STM4081'] == 'optimal':
                tmp_result = single_solution[0]
                tmp_result = floor(100*tmp_result['STM4081'])/100
                if tmp_result == tpiA_result:
                    print 'Passed MOMA single_deletion with tpiA & metN deletion in %1.4f seconds'%(time() - start_time)
            else:
                print 'failed single deletion'
            start_time - time()
            double_solution = double_deletion(wt_model, loci_list, loci_list, method='moma', solver=solver)
            tmp_result = floor(100*double_solution['data'][1,0])/100
            tmp_tpiA_result = floor(100*double_solution['data'][0,0])/100
            if tmp_result == tpiA_metN_result and tmp_tpiA_result == tpiA_result:
                print 'Passed MOMA double_deletion with tpiA & metN deletion in %1.4f seconds'%(time() - start_time)
            else:
                print 'failed double deletion'


#####
#turn into a double deletion unit test
if __name__ == '__main__':
    from cPickle import load
    from time import time
    from math import floor
    from numpy import array
    the_problem='return'
    element_type='gene'
    error_reporting=None
    from cobra.manipulation import initialize_growth_medium
    from cobra.test import salmonella_pickle
    with open(salmonella_pickle) as in_file:
        cobra_model = load(in_file)
    initialize_growth_medium(cobra_model, 'LB')
    the_names = ['tpiA', 'metN', 'atpA', 'eno']
    the_loci =  ['STM4081', 'STM0247', 'STM3867', 'STM2952']
    the_genes = tpiA, metN, atpA, eno = map(cobra_model.genes.get_by_id, the_loci)
    growth_dict = {'moma': {tpiA:1.61, metN:2.39, atpA:1.40, eno:0.33},
                   'fba':{tpiA:2.41, metN:2.43, atpA:1.87, eno:1.81}}
    growth_dict = {'moma':{'data': array([[1.61, 1.60, 0.95, 0.17],
                                          [1.60, 2.39, 1.40, 0.30],
                                          [0.95, 1.40, 1.40, 0.00],
                                          [0.17, 0.30, 0.00, 0.33]]),
                           'x': the_genes,
                           'y': the_genes},
                   'fba': {'data': array([[2.41, 2.38, 1.77, 1.81],
                                          [2.38, 2.43, 1.86, 1.79],
                                          [1.77, 1.86, 1.87, 1.32],
                                          [1.81, 1.79, 1.32, 1.81]]),
                           'x': the_genes,
                           'y': the_genes}
                   }
    solver_list = ['glpk',
                   'gurobi',
                   'cplex']
    try:
        import glpk
    except:
        solver_list.remove('glpk')
    try:
        from gurobipy import Model
    except:
        solver_list.remove('gurobi')
    try:
        from cplex import Cplex
    except:
        solver_list.remove('cplex')
 
    for solver in solver_list:
        print 'testing solver: ' + solver
        for method, the_growth_rates in growth_dict.items():
            if method == 'moma':
                print "MOMA isn't functional now, check back later"
                continue
            print '\twith method: ' + method 
            element_list_1 = the_growth_rates['x']
            element_list_2 = the_growth_rates['y']
            data = the_growth_rates['data']
            start_time = time()
            the_solution = double_deletion(cobra_model, element_list_1=element_list_1,
                                           element_list_2=element_list_2,
                                           method=method,  the_problem=the_problem,
                                           element_type=element_type,
                                           solver=solver,
                                           error_reporting=error_reporting)

            s_data = the_solution['data']
            print 'Double deletion simulation of %i genes ran in %1.3f seconds'%(s_data.size,
                                                                                 time()-start_time)
            s_x = the_solution['x']
            s_y = the_solution['y']
            for gene_x in element_list_1:
                for gene_y in element_list_2:
                    expected_value = data[element_list_1.index(gene_x),
                                          element_list_2.index(gene_y)]
                    simulated_value = floor(100*s_data[s_x.index(gene_x),
                                                       s_y.index(gene_y)])/100
                    if simulated_value != expected_value:
                        print '\t\tFAILED: %s/%s simulation (%1.3f) != expectation (%1.3f)'%(gene_x.name,
                                                                                             gene_y.name,
                                                                                             simulated_value,
                                                                                             expected_value)
                    else:
                        print '\t\tPASSED: %s/%s simulation (%1.3f) ~= expectation (%1.3f)'%(gene_x.name,
                                                                                             gene_y.name,
                                                                                             simulated_value,
                                                                                         expected_value)

            
