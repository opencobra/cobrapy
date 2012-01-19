#cobra.flux_analysis.essentiality.py
#runs flux variablity analysis on a Model object.
from warnings import warn
from math import floor,ceil
from numpy import vstack,zeros
from scipy import sparse
from copy import deepcopy
from cPickle import dump
from os import path, mkdir
try:
    #Allow for parallel simulations if ppmap is available
    from cobra.external import ppmap
    from double_deletion import double_deletion_parallel    
except:
    from double_deletion import double_deletion
    ppmap = False
from cobra.flux_analysis.moma import moma
from single_deletion import single_deletion
from cobra.manipulation import initialize_growth_medium
def assess_medium_component_essentiality(cobra_model, the_components=None,
                                         the_medium=None, solver='glpk',
                                         the_problem='return',
                                         the_condition=None, method='fba'):
    """Deterimes which components in an in silico medium are essential for growth in the
    context of the remaining components.

    cobra_model: A Model object.

    the_components: None or a list of exchange reactions that will be sequentially
    disabled.

    the_medium: Is None, a string, or a dictionary.  If a string then the
    initialize_growth_medium function expects that the_model has an
    attribute dictionary called media_compositions, which is a dictionary of
    dictionaries for various medium compositions.  Where a medium
    composition is a dictionary of exchange reaction ids for the medium
    components and the exchange fluxes for each medium component; note that
    these fluxes must be negative because they are being exchanged into the
    system.

    solver: 'glpk', 'gurobi', or 'cplex'

    the_problem: Is None, 'return', or an LP model object for the solver.
       
    the_condition: None or a String that provides a description of the medium
    simulation 

    """
    if method.lower() == 'moma':
        wt_model = cobra_model.copy()
    if isinstance(cobra_model, tuple):
        if len(cobra_model) == 3:
            the_condition = cobra_model[2]
        the_components = cobra_model[1]
        cobra_model = cobra_model[0]
    cobra_model = cobra_model.copy()
    if not the_components:
        if the_medium:
            if hasattr(the_medium, 'keys') or \
                   (hasattr(cobra_model,'media_compositions') and \
                    the_medium in cobra_model.media_compositions):
                initialize_growth_medium(cobra_model, the_medium)
                the_components = cobra_model.media_compositions[the_medium]
            else:
                raise Exception("%s is not a dict and not in the model's media list"%the_medium)
        else:
            raise Exception("You need to specify the_components or the_medium")
    essentiality_dict = {}
    for the_component in the_components:
        component_index = cobra_model.reactions.index(the_component)
        tmp_lb = float(cobra_model._lower_bounds[component_index])
        cobra_model.reactions[component_index].lower_bound = cobra_model._lower_bounds[component_index] = 0
        if method.lower() == 'fba':
            cobra_model.optimize(solver=solver, the_problem=the_problem)
            objective_value = cobra_model.solution.f
        elif method.lower() == 'moma':
           objective_value = moma(wt_model, cobra_model, solver=solver)['objective_value'] 
        essentiality_dict[the_component] = objective_value
        cobra_model.reactions[component_index].lower_bound = cobra_model._lower_bounds[component_index] = tmp_lb
    if the_condition:
        essentiality_dict['the_condition'] = the_condition
    return(essentiality_dict)

def deletion_analysis(cobra_model, the_medium=None, deletion_type='single',
                      work_directory=None, growth_cutoff=0.001,
                      the_problem='return', n_processes=6, element_type='gene',
                      solver='glpk', error_reporting=None, method='fba', element_list=None):
    """Performs single and/or double deletion analysis on all the genes in the model.  Provides
    an interface to parallelize the deletion studies.

    cobra_model: A Model object.

    the_medium: Is None, a string, or a dictionary.  If a string then the
    initialize_growth_medium function expects that cobra_model has an
    attribute dictionary called media_compositions, which is a dictionary of
    dictionaries for various medium compositions.  Where a medium
    composition is a dictionary of exchange reaction ids for the medium
    components and the exchange fluxes for each medium component; note that
    these fluxes must be negative because they are being exchanged into the
    system.

    deletion_type: 'single', 'double', or 'double-only'

    work_directory: None or String indicating where to save the output from the simulations.

    growth_cutoff: Float.  Indiates the minimum growth rate that is considered viable.

    the_problem: Is None, 'return', or an LP model object for the solver.
       
    element_type: 'gene' or 'reaction'

    solver: 'glpk', 'gurobi', or 'cplex'

    n_processes: number of parallel processes to break the double deletion simulations
    into.

    error_reporting: None or True

    element_list: None or a list of genes to delete from the model.

    Returns: Nothing.  However, the script will add attributes single_deletion_* and
    double_deletion_* to cobra_model containing the simulation results.

    """
    if element_type == 'reaction':
        warn("deletion_analysis is not perfect for element_type = 'reaction'")
    #When using ppmap, it's easier to feed in the parameters as a list,
    #if the defaults need to be changed
    if isinstance(cobra_model, list):
        tmp_model = cobra_model
        cobra_model = tmp_model[0]
        if len(tmp_model) > 1:
            the_medium = tmp_model[1]
            if len(tmp_model) > 2:
                deletion_type = tmp_model[2]
                if len(tmp_model) > 3:
                    work_directory = tmp_model[3]
                    if len(tmp_model) > 4:
                        growth_cutoff = tmp_model[4]
    if the_medium is not None:
        initialize_growth_medium(cobra_model, the_medium)
    the_problem=cobra_model.optimize(the_problem=the_problem, solver=solver)
    #Store the basal model for the simulations
    if element_list is None:
        element_list = getattr(cobra_model, element_type + 's')
    if deletion_type != 'double_only':
        cobra_model.single_deletion_growth_wt = cobra_model.solution.f
        growth_rate_dict, growth_solution_status_dict, problem_dict = single_deletion(deepcopy(cobra_model),
                                                                                      element_list=element_list,
                                                                                      the_problem=the_problem,
                                                                                      element_type=element_type,
                                                                                      solver=solver,
                                                                                      error_reporting=error_reporting,
                                                                                      method=method)
        del problem_dict
        cobra_model.single_deletion_growth_dict = growth_rate_dict
        cobra_model.single_deletion_solution_status_dict = growth_solution_status_dict
        setattr(cobra_model, 'single_deletion_%ss'%element_type, deepcopy(growth_rate_dict.keys()))
        cobra_model.single_deletion_lethal = [x for x in growth_rate_dict.keys()
                                            if growth_rate_dict[x] <  growth_cutoff]

        cobra_model.single_deletion_growth_medium = the_medium 
        cobra_model.single_deletion_nonlethal =  list(set(growth_rate_dict.keys()).difference(cobra_model.single_deletion_lethal))
        if work_directory is not None:
            if not path.lexists(work_directory):
                mkdir(work_directory)
            with open(work_directory + the_medium + '_single_' + cobra_model.description, 'w') as out_file:
                dump(cobra_model, out_file)

    if deletion_type == 'double' or deletion_type == 'double_only':
        #It appears that the glpk interface no longer works will with sending
        #a glpk.LPX object through ppmap, so just set the basis to return
        if the_problem:
            the_problem='return'
        cobra_model.double_deletion_growth_medium = the_medium 
        cobra_model.double_deletion_growth_wt = cobra_model.solution.f
        if not ppmap:
            if n_processes > 0:
                print "Couldn't import ppmap from cobra.external is parallel python installed?"
                return

        else:
            cobra_model = double_deletion_parallel(deepcopy(cobra_model),
                                                 genes_of_interest=element_list,
                                                 the_problem=the_problem,
                                                 n_processes=n_processes,
                                                 element_type=element_type,
                                                 solver=solver,
                                                 error_reporting=error_reporting,
                                                 method=method)
        #This indicates the genes that were run through double deletion but
        #the x and y lists specify the order
        setattr(cobra_model, 'double_deletion_%ss'%element_type, deepcopy(cobra_model.genes))
        if work_directory is not None:
            with open(work_directory + the_medium + '_double_' + cobra_model.description, 'w') as out_file:
                dump(cobra_model, out_file)
    return cobra_model 
