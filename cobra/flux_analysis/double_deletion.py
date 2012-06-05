#cobra.flux_analysis.double_deletion.py
#runs flux variablity analysis on a Model object.
from pdb import set_trace
from math import floor,ceil
from numpy import vstack,zeros, nan
from scipy import sparse
from copy import deepcopy
from cobra.flux_analysis.moma import moma
from warnings import warn
try:
    from cobra.external.ppmap import ppmap
except:
    ppmap = False
    print 'Parallel mode is not available.  Is Parallel Python installed?'
from cPickle import dump
#from double_deletion import double_deletion
from cobra.manipulation import initialize_growth_medium
from cobra.manipulation import delete_model_genes, undelete_model_genes
def double_deletion_parallel(cobra_model, n_processes=4,
                             genes_of_interest=None, method='fba', the_medium=None,
                             the_problem='return', element_type='gene',
                             solver='glpk',
                             error_reporting=None):
    """Wrapper for double_gene_deletion_parallel and the currently
    unimplemented double_reaction_deletion_parallel functions

    cobra_model: a cobra.Model object

    n_processes: is the number of parallel processes to start

    genes_of_interest: Is None, a list of genes, or a list of two lists of
    genes.  If None then double_deletion is run on all genes in
    the_model.genes.  If a list of genes then double_deletion is run for all
    combinations of genes in double_deletion.  If a list of of two lists of
    genes then double_deletion is run for each member of one list vs. each
    member of the second list.

    method: 'fba' or 'moma' to run flux balance analysis or minimization
    of metabolic adjustments.
    
    the_medium: Is None, a string, or a dictionary.  If a string then the
    initialize_growth_medium function expects that the_model has an
    attribute dictionary called media_compositions, which is a dictionary of
    dictionaries for various medium compositions.  Where a medium
    composition is a dictionary of exchange reaction ids for the medium
    components and the exchange fluxes for each medium component; note that
    these fluxes must be negative because they are being exchanged into the
    system.

    the_problem: Is None or 'reuse'

    element_type: 'gene' or 'reaction'

    solver: 'glpk', 'gurobi', or 'cplex'.

    error_reporting: None or True.

    Returns a dictionary of the genes in the x dimension (x), the y
    dimension (y), and the growth simulation data (data).

    """
    if not ppmap:
        if not ppmap:
           print  'Parallel mode not available is Parallel Python installed?'
           return
    if the_problem:
        #The solver model objects are not thread safe so change the_problem
        #to 'return' if one is provided.
        the_problem='return'
    if element_type == 'gene':
        return double_gene_deletion_parallel(cobra_model, n_processes=n_processes,
                                             genes_of_interest=genes_of_interest,
                                             method=method, the_medium=the_medium,
                                             the_problem=the_problem, solver=solver,
                                             error_reporting=error_reporting)
    else:
        raise Exception("Double deletion not yet implemented for element_type = %s"%element_type)


def double_deletion(cobra_model, element_list_1=None, element_list_2=None,
                    method='fba', single_deletion_growth_dict=None, the_problem='return',
                    element_type='gene', solver='glpk', error_reporting=None):
    """Wrapper for double_gene_deletion and the currently unimplemented
    double_reaction_deletion functions

    cobra_model: a cobra.Model object

    element_list_1: Is None or a list of elements (genes or reactions)
    
    element_list_2: Is None or a list of elements (genes or reactions)

    method: 'fba' or 'moma' to run flux balance analysis or minimization
    of metabolic adjustments.

    single_deletion_growth_dict: A dictionary that provides the growth
    rate information for single gene knock outs.  This can speed up
    simulations because nonviable single deletion strains imply that all
    double deletion strains will also be nonviable.

    the_problem: Is None or 'reuse'

    element_type: 'gene' or 'reaction'

    solver: 'glpk', 'gurobi', or 'cplex'.

    error_reporting: None or True

    Returns a dictionary of the elements in the x dimension (x), the y
    dimension (y), and the growth simulation data (data).

    """
    if element_type == 'gene':
        return double_gene_deletion(cobra_model, gene_list_1=element_list_1,
                                    gene_list_2=element_list_2, method=method,
                                    single_deletion_growth_dict=single_deletion_growth_dict,
                                    the_problem=the_problem, solver=solver,
                                    error_reporting=error_reporting)
    else:
        raise Exception("Double deletion not yet implemented for element_type = %s"%element_type)

def double_gene_deletion(cobra_model, gene_list_1=None, gene_list_2=None,
                         method='fba', single_deletion_growth_dict=None,
                         the_problem='return', solver='glpk',
                         error_reporting=None):
    """This will disable reactions for all gene pairs from gene_list_1 and
    gene_list_2 and then run simulations to optimize for the objective
    function.  The contribution of each reaction to the objective function
    is indicated in cobra_model._objective_coefficients vector.

    cobra_model: a cobra.Model object

    gene_list_1: Is None or a list of genes.  If None then both gene_list_1
    and gene_list_2 are assumed to correspond to cobra_model.genes.
    
    gene_list_2: Is None or a list of genes.  If None then gene_list_2 is
    assumed to correspond to gene_list_1.

    method: 'fba' or 'moma' to run flux balance analysis or minimization
    of metabolic adjustments.
    
    single_deletion_growth_dict: A dictionary that provides the growth
    rate information for single gene knock outs.  This can speed up
    simulations because nonviable single deletion strains imply that all
    double deletion strains will also be nonviable.

    the_problem: Is None, 'return', or an LP model object for the solver.

    solver: 'glpk', 'gurobi', or 'cplex'.

    error_reporting: None or True

    Returns a dictionary of the genes in the x dimension (x), the y
    dimension (y), and the growth simulation data (data).
    
    """
    from warnings import warn

    warn("Need to use more of cobra.Gene and less of strs in gene_list")
    #BUG: Since this might be called from ppmap, the modules need to
    #be imported.  Modify ppmap to take depfuncs
    from numpy import zeros, nan
    from cobra.flux_analysis.single_deletion import single_deletion
    from cobra.manipulation import initialize_growth_medium
    from cobra.manipulation import delete_model_genes, undelete_model_genes
    ##TODO: Use keywords instead
    if isinstance(cobra_model, dict):
        tmp_dict = cobra_model
        cobra_model = tmp_dict['cobra_model']
        if 'gene_list_1' in tmp_dict:
            gene_list_1 = tmp_dict['gene_list_1']
        if 'gene_list_2' in tmp_dict:
            gene_list_2 = tmp_dict['gene_list_2']
        if 'method' in tmp_dict:
            method = tmp_dict['method']
        if 'the_problem' in tmp_dict:
            the_problem = tmp_dict['the_problem']
        if 'single_deletion_growth_dict' in tmp_dict:
            single_deletion_growth_dict = tmp_dict['single_deletion_growth_dict']
        if 'solver' in tmp_dict:
            solver = tmp_dict['solver']
        if 'error_reporting' in tmp_dict:
            error_reporting = tmp_dict['error_reporting']
    else:
        cobra_model = cobra_model
    #this is a slow way to revert models.
    wt_model = cobra_model  #NOTE: It may no longer be necessary to use a wt_model
    #due to undelete_model_genes
    if gene_list_1 is None:
        gene_list_1 = cobra_model.genes
    elif not hasattr(gene_list_1[0], 'id'):
        gene_list_1 = map(cobra_model.genes.get_by_id, gene_list_1)
    #Get default values to use if the deletions do not alter any reactions
    the_problem = cobra_model.optimize(the_problem=the_problem, solver=solver)
    basal_f = cobra_model.solution.f
    if method.lower() == 'moma':
        wt_model = cobra_model.copy()
        the_problem = 'return'
        combined_model = None
    single_gene_set = set(gene_list_1)
    if gene_list_2:
        if not hasattr(gene_list_2[0], 'id'):
            gene_list_2 = map(cobra_model.genes.get_by_id, gene_list_2)
        single_gene_set.update(gene_list_2)
    #Run the single deletion analysis to account for double deletions that
    #target the same gene and lethal deletions.  We assume that there
    #aren't synthetic rescues.
    single_deletion_growth_dict = single_deletion(cobra_model,
                                                  list(single_gene_set),
                                                  method=method,
                                                  the_problem=the_problem,
                                                  solver=solver,
                                                  error_reporting=error_reporting)[0]
    if gene_list_2 is None or gene_list_1 == gene_list_2:
        deletion_array = zeros([len(gene_list_1), len(gene_list_1)]) 
        ##TODO: Speed up this triangular process
        #For the case where the contents of the lists are the same cut the work in half.
        #There might be a faster way to do this by using a triangular array function
        #in numpy
        #Populate the diagonal from the single deletion lists
        for i, the_gene in enumerate(gene_list_1):
            deletion_array[i, i] = single_deletion_growth_dict[the_gene]
        for i in range(len(gene_list_1)-1):
            gene_1 = gene_list_1[i]
            #TODO: Since there cannot be synthetic rescues we can assume
            #that the whole row for a lethal deletion
            #will be equal to that deletion.
            if single_deletion_growth_dict[gene_1] <= 0:
                for j in range(i+1, len(gene_list_1)):
                     deletion_array[j, i] = deletion_array[i, j] = single_deletion_growth_dict[gene_2]
            else:
                for j in range(i+1, len(gene_list_1)):
                    if single_deletion_growth_dict[gene_1] <= 0:
                        tmp_solution = single_deletion_growth_dict[gene_1]
                    else:
                        gene_2 = gene_list_1[j]
                        delete_model_genes(cobra_model, [gene_1, gene_2])
                        if cobra_model._trimmed:
                            if method.lower() == 'fba':
                                #Assumes that the majority of perturbations don't change
                                #reactions which is probably false
                                cobra_model.optimize(the_problem = the_problem,
                                                     solver=solver, error_reporting=error_reporting)
                                the_status = cobra_model.solution.status
                                tmp_solution = cobra_model.solution.f
                            elif method.lower() == 'moma':
                                try:
                                    moma_solution = moma(wt_model, cobra_model,
                                                         combined_model=combined_model,
                                                         solver=solver, the_problem=the_problem)
                                    tmp_solution = float(moma_solution.pop('objective_value'))
                                    the_problem = moma_solution.pop('the_problem')
                                    the_status = moma_solution.pop('status')
                                    combined_model = moma_solution.pop('combined_model')
                                    del moma_solution
                                except:
                                    tmp_solution = nan
                                    the_status = 'failed'
                            if the_status not in ['opt', 'optimal']  and \
                                   error_reporting:
                                print '%s / %s: %s status: %s'%(gene_1, gene_2, solver,
                                                                the_status)   
                            #Reset the model to orginial form.
                            undelete_model_genes(cobra_model)
                        else:
                            tmp_solution = basal_f
                    deletion_array[j, i] = deletion_array[i, j] = tmp_solution

    else:
        deletion_array = zeros([len(gene_list_1), len(gene_list_2)])
        #Now deal with the case where the gene lists are different
        for i, gene_1 in enumerate(gene_list_1):
            if single_deletion_growth_dict[gene_1] <= 0:
                for j in range(len(gene_list_2)):
                    deletion_array[i, j] = 0.
            else:
                for j, gene_2 in enumerate(gene_list_2):
                    #Assume no such thing as a synthetic rescue
                    if single_deletion_growth_dict[gene_2] <= 0:
                        tmp_solution = single_deletion_growth_dict[gene_2]
                    else:
                        delete_model_genes(cobra_model, [gene_1, gene_2])
                        if cobra_model._trimmed:
                            if method.lower() == 'fba':
                                cobra_model.optimize(the_problem=the_problem,
                                                     solver=solver, error_reporting=error_reporting)
                                tmp_solution = cobra_model.solution.f
                                the_status = cobra_model.solution.status
                            elif method.lower() == 'moma':
                                try:
                                    moma_solution = moma(wt_model, cobra_model,
                                                         combined_model=combined_model,
                                                         solver=solver, the_problem=the_problem)
                                    tmp_solution = float(moma_solution.pop('objective_value'))
                                    the_problem = moma_solution.pop('the_problem')
                                    the_status = moma_solution.pop('status')
                                    combined_model = moma_solution.pop('combined_model')
                                    del moma_solution
                                except:
                                    tmp_solution = nan
                                    the_status = 'failed'
                            if the_status not in ['opt', 'optimal']  and \
                                   error_reporting:
                                print '%s / %s: %s status: %s'%(repr(gene_1), repr(gene_2), solver,
                                                            cobra_model.solution.status)
                            #Reset the model to wt form
                            undelete_model_genes(cobra_model)
                        else:
                            tmp_solution = basal_f
                    deletion_array[i, j] = tmp_solution

    return({'x': gene_list_1, 'y': gene_list_2, 'data': deletion_array})


def double_gene_deletion_parallel(cobra_model, n_processes=4,
                                  genes_of_interest=None, method = 'fba', the_medium=None,
                                  the_problem='return', solver='glpk',
                                  error_reporting=None):
    """Provides a wrapper to run the double_deletion function on
    multicore systems.

    cobra_model: a Model object

    n_processes: is the number of parallel processes to start

    genes_of_interest: Is None, a list of genes, or a list of two lists of
    genes.  If None then double_deletion is run on all genes in
    cobra_model.genes.  If a list of genes then double_deletion is run for all
    combinations of genes in double_deletion.  If a list of of two lists of
    genes then double_deletion is run for each member of one list vs. each
    member of the second list.

    method: 'fba' or 'moma' to run flux balance analysis or minimization
    of metabolic adjustments.

    the_medium: Is None, a string, or a dictionary.  If a string then the
    initialize_growth_medium function expects that cobra_model has an
    attribute dictionary called media_compositions, which is a dictionary of
    dictionaries for various medium compositions.  Where a medium
    composition is a dictionary of exchange reaction ids for the medium
    components and the exchange fluxes for each medium component; note that
    these fluxes must be negative because they are being exchanged into the
    system.

    the_problem: Is None or 'reuse'

    solver: 'glpk', 'gurobi', or 'cplex'.

    error_reporting: None or True

    Adds the following attributes to the cobra_model:
          double_deletion_growth_rate: A numpy array of the simulation results
          double_deletion_genes_x: A list of the genes for the x dimension of
          double_deletion_growth_rate.
          double_deletion_genes_y: A list of the genes for the y dimension of
          double_deletion_growth_rate.

    """
    from warnings import warn

    warn("Need to use more of cobra.Gene and less of strs in gene_list")
    if not ppmap:
        print  'Parallel mode not available is Parallel Python installed'
        return
    if the_problem:
        the_problem='return'
    if the_medium:
        initialize_growth_medium(cobra_model, the_medium)
    if not genes_of_interest:
        #If no genes_of_interest are specified then assume we want to
        #compare all genetic interactions in the network
        all_genes = [x.id for x in cobra_model.genes]
        second_gene_list = all_genes
    elif isinstance(genes_of_interest[0], str):
        #If genes_of_interest is a list then assume the list be scanned
        #for interactions with all genes in the network
        all_genes = genes_of_interest
        second_gene_list = all_genes
    elif hasattr(genes_of_interest[0], 'id'):
        #Make sure we're dealing with strings instead of objects because we
        #haven't audited this for thread safety
        second_gene_list = all_genes = [x.id for x in genes_of_interest]
    elif hasattr(genes_of_interest[0], '__iter__'):
        all_genes = genes_of_interest[0]
        if len(genes_of_interest) == 2:
            second_gene_list = genes_of_interest[1] 
        else:
            second_gene_list = all_genes
    #Get basic numbers to guide how the problem should be divided for parallel execution.
    total_gene_count = len(all_genes)
    division_count = total_gene_count / n_processes
    the_rows = []
    #In the future implement passing the single del list to each spawned double
    #del process to save some time.
    # single_deletion_growth_dict, tmp_stat = tools.single_deletion( deepcopy( cobra_model),
    #gene_list = list( set( all_genes ).union( second_gene_list ) ), the_problem = the_problem )
    for i in range(n_processes-1):
        the_rows.append({'cobra_model': cobra_model.copy(), 'method': method,
                         'gene_list_1': deepcopy(all_genes[i*division_count:division_count*(i+1)]),
                         'gene_list_2': deepcopy(second_gene_list), 'the_problem': the_problem,
                         'solver': solver,
                         'error_reporting': error_reporting})
    the_rows.append({'cobra_model': cobra_model.copy(), 'method': method,
                     'gene_list_1': deepcopy(all_genes[(n_processes-1)*division_count:]),
                     'gene_list_2': deepcopy(second_gene_list), 'the_problem': the_problem,
                     'solver': solver,
                     'error_reporting': error_reporting})

    tmp_pp = list(ppmap(n_processes, double_gene_deletion, the_rows))
    gene_list_x = tmp_pp[0]['x']
    gene_list_y = tmp_pp[0]['y']
    double_deletion_data = tmp_pp[0]['data']
    for the_result in tmp_pp[1:]:
        gene_list_x += the_result['x']
        double_deletion_data = vstack((double_deletion_data, the_result['data']))

    cobra_model.double_deletion_growth_rate = double_deletion_data 
    cobra_model.double_deletion_genes_x = gene_list_x
    cobra_model.double_deletion_genes_y = gene_list_y
    return cobra_model


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
