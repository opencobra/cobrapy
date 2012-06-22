#cobra.flux_analysis.moma.py: Runs the minimization of metabolic
#adjustment method described in Segre et al 2002 PNAS 99(23): 15112-7
from copy import deepcopy
from time import time
from math import ceil, floor
from numpy import array, zeros, ones, hstack, vstack, matrix, sum
from scipy.sparse import eye, lil_matrix, dok_matrix
from scipy.sparse import hstack as s_hstack
from scipy.sparse import vstack as s_vstack
from ..core import Reaction, Metabolite
#from cobra.core.Metabolite import Metabolite
from ..manipulation import initialize_growth_medium, delete_model_genes
from warnings import warn
#TODO: Add in an option for using matrices instead of objects because it
#appears that there might be a performance penalty (especially for repetitions)
#when using objects.
#
#
#Using solver: cplex with lp_method: 1
#
#older optimized matrix method:
#    Passed MOMA with tpiA deletion in 5.0713 seconds
#    Passed MOMA reusing Model with tpiA deletion in 4.0048 seconds
#    Passed MOMA reusing Model and model with tpiA deletion in 0.6579 seconds
#    Passed MOMA single_deletion with tpiA & metN deletion in 14.6137 seconds
#    Passed MOMA double_deletion with tpiA & metN deletion in 41.2922 seconds
#
#new unoptimized object method:
#    add time 7.898585
#    Took 8.954579 seconds to construct combined model
#    Took 0.017424 seconds to update combined model
#    Took 1.043317 seconds to solve problem
#    Took 0.002228 seconds to assemble solution
#    Passed MOMA with tpiA deletion in 10.1922 seconds
#    Passed MOMA single_deletion with tpiA & metN deletion in 30.0383 seconds
#    Passed MOMA double_deletion with tpiA & metN deletion in 69.8822 seconds
#
#  The major penalties are related to adding two models (cobra.core.Model.__add__)


def moma(wt_model, mutant_model, objective_sense='maximize', solver='gurobi',
         tolerance_optimality=1e-8, tolerance_feasibility=1e-8,
         minimize_norm=False, the_problem='return', lp_method=0,
         combined_model=None, norm_type='euclidean', print_time=False):
    """Runs the minimization of metabolic adjustment method described in
    Segre et al 2002 PNAS 99(23): 15112-7.

    wt_model: A cobra.Model object

    mutant_model: A cobra.Model object with different reaction bounds vs wt_model.
    To simulate deletions

    objective_sense: 'maximize' or 'minimize'

    solver: 'gurobi', 'cplex', or 'glpk'.  Note: glpk cannot be used with
    norm_type 'euclidean'

    tolerance_optimality: Solver tolerance for optimality.

    tolerance_feasibility: Solver tolerance for feasibility.

    the_problem: None or a problem object for the specific solver that can be
    used to hot start the next solution.

    lp_method: The method to use for solving the problem.  Depends on the solver.  See
    the cobra.flux_analysis.solvers.py file for more info.
        For norm_type == 'euclidean':
            the primal simplex works best for the test model (gurobi: lp_method=0, cplex: lp_method=1)
    
    combined_model: an output from moma that represents the combined optimization
    to be solved.  when this is not none.  only assume that bounds have changed
    for the mutant or wild-type.  This saves 0.2 seconds in stacking matrices.

    """
    warn('MOMA is currently non-functional.  check back later')
    if solver.lower() == 'cplex' and lp_method == 0:
        #print 'for moma, solver method 0 is very slow for cplex. changing to method 1'
        lp_method = 1
    if solver.lower() == 'glpk' and norm_type == 'euclidean':
        print "GLPK can't solve quadratic problems like MOMA.  Switching to linear MOMA"

    if norm_type == 'euclidean':
        #Reusing the basis can get the solver stuck.
        reuse_basis = False
    if combined_model and combined_model.norm_type != norm_type:
        print 'Cannot use combined_model.norm_type = %s with user-specified norm type'%(combined_model.norm_type,
                                                                                        norm_type)
        print 'Defaulting to user-specified norm_type'
        combined_model = None
    #Add a prefix in front of the mutant_model metabolites and reactions to prevent
    #name collisions in DictList
    for the_dict_list in [mutant_model.metabolites,
                          mutant_model.reactions]:
        [setattr(x, 'id', 'mutant_%s'%x.id)
         for x in the_dict_list]
        the_dict_list._generate_index() #Update the DictList.dicts


    wt_model.optimize(solver=solver)
    wt_solution = deepcopy(wt_model.solution)
    if objective_sense == 'maximize':
        wt_optimal = floor(wt_solution.f/tolerance_optimality)*tolerance_optimality
    else:
        wt_optimal = ceil(wt_solution.f/tolerance_optimality)*tolerance_optimality
    if norm_type == 'euclidean':
        quadratic_component = eye(wt_solution.x.shape[0],wt_solution.x.shape[0])
    elif norm_type == 'linear':
        raise Exception('linear MOMA is not currently implmented')
        quadratic_component = None
    if minimize_norm:
        raise Exception('minimize_norm is not currently implemented')
        #just worry about the flux distribution and not the objective from the wt
        combined_model = mutant_model.copy()
        #implement this: combined_model.reactions[:].objective_coefficients = -wt_solution.x_dict
    else:
        #Construct a problem that attempts to maximize the objective in the WT model while
        #solving the quadratic problem.  This new problem is constructed to try to find
        #a solution for the WT model that lies close to the mutant model.  There are
        #often multiple equivalent solutions with M matrices and the one returned
        #by a simple cobra_model.optimize call may be too far from the mutant.
        #This only needs to be adjusted if we update mutant_model._S after deleting reactions

        if print_time:
            start_time = time()
        number_of_reactions = len(mutant_model.reactions)
        if norm_type == 'euclidean':
            reaction_coefficient = 1
        elif norm_type == 'linear':
            reaction_coefficient = 2
        if not combined_model:
            #Collect the set of wt reactions contributing to the objective.
            objective_reaction_coefficient_dict = dict([(x.id, x.objective_coefficient)
                                                        for x in wt_model.reactions
                                                        if x.objective_coefficient])
            #This does a deepcopy of both models which might result in a huge overhead.
            #Update cobra.core.Model to improve performance.
            combined_model = wt_model + mutant_model
            if print_time:
                print 'add time %f'%(time()-start_time)
            [setattr(x, 'objective_coefficient', 0.)
             for x in combined_model.reactions]
            #Add in the difference reactions.  The mutant reactions and metabolites are already added.
            #This must be a list to maintain the correct order when adding the difference_metabolites

            difference_reactions = [Reaction('difference_%i'%i)
                                        for i in range(reaction_coefficient*number_of_reactions)]
            [setattr(x, 'lower_bound', -1000)
             for x in difference_reactions]
            combined_model.add_reactions(difference_reactions)
            index_to_reaction = combined_model.reactions
            id_to_reaction = combined_model.reactions._object_dict
            #This is slow
            #Add in difference metabolites
            difference_metabolite_dict = dict([(i, Metabolite('difference_%i'%i))
                                           for i in xrange(number_of_reactions)])
            combined_model.add_metabolites(difference_metabolite_dict.values())
            for i, tmp_metabolite in difference_metabolite_dict.iteritems():
                if norm_type == 'linear':
                    tmp_metabolite._constraint_sense = 'G'
                index_to_reaction[i].add_metabolites({tmp_metabolite: -1.},
                                                     add_to_container_model=False)
                index_to_reaction[i+number_of_reactions].add_metabolites({tmp_metabolite: 1.}, add_to_container_model=False)
                index_to_reaction[i+2*number_of_reactions].add_metabolites({tmp_metabolite: 1.}, add_to_container_model=False)

            #Add in the virtual objective metabolite
            objective_metabolite = Metabolite('wt_optimal')
            objective_metabolite._bound = wt_optimal
            if objective_sense == 'maximize':
                objective_metabolite._constraint_sense = 'G'
            else:
                objective_metabolite._constraint_sense = 'L'
            #TODO: this couples the wt_model objective reaction to the virtual metabolite
            #Currently, assumes a single objective reaction; however, this may be extended
            [id_to_reaction[k].add_metabolites({objective_metabolite: v})
             for k, v in objective_reaction_coefficient_dict.items()]

            if print_time:
                print 'Took %f seconds to construct combined model'%(time()-start_time)
                start_time = time()



        if norm_type == 'euclidean':
            quadratic_component = s_vstack((lil_matrix((2*number_of_reactions, 3*number_of_reactions)),
                                            s_hstack((lil_matrix((number_of_reactions, 2*number_of_reactions)),
                                                      quadratic_component))))
        elif norm_type == 'linear':
            quadratic_component = None

    combined_model.norm_type = norm_type
    cobra_model = combined_model

    if print_time:
        print 'Took %f seconds to update combined model'%(time()-start_time)
        start_time = time()
    the_result = combined_model.optimize(objective_sense='minimize',
                                         quadratic_component=quadratic_component,
                                         solver=solver,
                                         tolerance_optimality=tolerance_optimality,
                                         tolerance_feasibility=tolerance_feasibility,
                                         lp_method=lp_method, reuse_basis=reuse_basis)
    the_problem = the_result
    the_solution = combined_model.solution

    if print_time:
        print 'Took %f seconds to solve problem'%(time()-start_time)
        start_time = time()
    mutant_dict = {}
    x_vector = the_solution.x
    if hasattr(x_vector, 'flatten'):
        x_vector = x_vector.flatten()
    mutant_dict['x'] = mutant_fluxes = array(x_vector[1*number_of_reactions:2*number_of_reactions])
    #Need to use the new solution as there are multiple ways to achieve an optimal solution in
    #simulations with M matrices.
    wt_model.solution.x = array(x_vector[:number_of_reactions])
    mutant_dict['objective_value'] = mutant_f = float(matrix(mutant_fluxes)*matrix([x.objective_coefficient
                                                                                    for x in mutant_model.reactions]).T)
    mutant_dict['status'] = the_solution.status
    mutant_dict['flux_difference'] = flux_difference = sum((wt_model.solution.x - mutant_fluxes)**2)
    mutant_dict['the_problem'] = the_problem
    mutant_dict['combined_model'] = combined_model
    if print_time:
        print 'Took %f seconds to assemble solution'%(time()-start_time)
    
    del wt_model, mutant_model, quadratic_component, the_solution
    return(mutant_dict)


