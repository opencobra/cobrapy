# -*- coding: utf-8 -*-
# PLEASE NOTE THAT JYTHON SUPPORT (and this jython-only-solver) is deprecated
#This script provides wrappers for libglpk-java 1.0.22 and pyglpk 0.3
from __future__ import absolute_import, print_function

from copy import deepcopy
from os import name
from time import time
from warnings import warn

from six import iteritems

from org.gnu.glpk import GLPK, GLPKConstants, glp_iocp, glp_smcp

from ..core.solution import Solution
###solver specific parameters
from .parameters import (
    default_objective_sense, objective_senses, parameter_defaults,
    parameter_mappings, sense_dict, status_dict, variable_kind_dict)

solver_name = 'glpk'
sense_dict = eval(sense_dict[solver_name])
#Functions that are different for java implementation of a solver

if name != "java":
    raise Exception("jython only")

warn("cobra.solvers.glpk_solver isn't mature.  consider using gurobi or cplex")
variable_kind_dict = eval(variable_kind_dict['%s_%s'%(solver_name,
                                            __name)])
status_dict = eval(status_dict['%s_%s'%(solver_name,
                                            __name)])
objective_senses = objective_senses['%s_%s'%(solver_name,
                                            __name)]
parameter_mappings = parameter_mappings['%s_%s'%(solver_name,
                                                 __name)]
parameter_defaults = parameter_defaults['%s_%s'%(solver_name,
                                                 __name)]

class Problem():
    """Create a more pythonesqe class to wrap the key
    features of the libglpk-java functions.
    
    """
    def __init__(self):
        """the attributes g, lp, mip should be made private
        """
        self._g = GLPK
        self._lp= GLPK.glp_create_prob()
        self._simplex_parameters = glp_smcp()
        self._mip_parameters = None
        self._g.glp_init_smcp(self._simplex_parameters)
        self.status = self.objective_value = None
        self._mip = False
    def set_name(self, name=''):
        self._g.glp_set_prob_name(self._lp, name)

    def solve(self):
        try:
            self._g.glp_simplex(self._lp,
                               self._simplex_parameters)
            if self._mip:
                #perform the MIP
                setattr(self._mip_parameters, 'msg_lev',
                         self._simplex_parameters.msg_lev)
                self._g.glp_intopt(self._lp, self._mip_parameters)
            self.status = self.get_status()
            self.objective_value = self.get_objective_value()
        except:
            self.status = 'failed'
        return self.status

    def get_status(self):
        if self._mip:
            status = self._g.glp_mip_status(self._lp)
        else:
            status = self._g.glp_get_status(self._lp)
        return status_dict[status]
    
    def set_objective_sense(self, parameter_value='maximize'):
        self._g.glp_set_obj_dir(self._lp,
                               eval(objective_senses[parameter_value]))

    def set_parameter(self, parameter_name, parameter_value, warning=False):
        if parameter_name == 'objective_sense':
            self.set_objective_sense(parameter_value)
        else:
            if parameter_name == 'meth' and parameter_value not in [1,2,3]:
                parameter_value = 1
            try:
                setattr(self._simplex_parameters, parameter_name,
                        parameter_value)
            except Exception as e1:
                try:
                    setattr(self._mip_parameters, parameter_name,
                            parameter_value)
                except Exception as e2:
                    if warning:
                        print("Could not set simplex parameter " +\
                              "{:s}: {:s}".format(parameter_name, repr(e1)))
                        
                        if self._mip_parameters is not None:
                            print("Could not set mip parameter " +\
                                  "{:s}: {:s}".format(parameter_name, repr(e2)))
    def get_objective_value(self):
        if self._mip:
            tmp_value = self._g.glp_mip_obj_val(self._lp)
        else:
            tmp_value = self._g.glp_get_obj_val(self._lp)
        return tmp_value

    def create_problem(self, cobra_model):
        g = self._g
        lp = self._lp
        number_of_reactions = len(cobra_model.reactions)
        number_of_metabolites = len(cobra_model.metabolites)
        g.glp_add_cols(lp, number_of_reactions)
        reaction_to_index = {}
        objective_dict = {}
        #Add in the variables
        tmp_kinds = []
        for i, the_reaction in enumerate(cobra_model.reactions):
            i_offset = i + 1
            reaction_to_index[the_reaction] = i_offset
            if the_reaction.objective_coefficient != 0:
                objective_dict[i_offset] = the_reaction.objective_coefficient
            g.glp_set_col_name(lp, i_offset, the_reaction.id)
            tmp_kinds.append(the_reaction.variable_kind)
            the_kind = variable_kind_dict[the_reaction.variable_kind]
            lower_bound = the_reaction.lower_bound
            upper_bound = the_reaction.upper_bound
            #Note. It is possible to have unbounded or one-bound variables
            if lower_bound == upper_bound:
                bound_kind = GLPKConstants.GLP_FX
            else:
                bound_kind = GLPKConstants.GLP_DB
            g.glp_set_col_kind(lp, i_offset, the_kind)
            g.glp_set_col_bnds(lp, i_offset,
                               bound_kind, the_reaction.lower_bound,
                               the_reaction.upper_bound)
        tmp_kinds = set(tmp_kinds)
        if 'integer' in tmp_kinds or 'binary' in tmp_kinds:
            self._mip = True
            self._mip_parameters = glp_iocp()
            g.glp_init_iocp(self._mip_parameters)
        #create constraints
        g.glp_add_rows(lp, number_of_metabolites)
        row_indices = []
        column_indices = []
        constraint_values = []
        for i, the_metabolite in enumerate(cobra_model.metabolites):
            i_offset = i + 1
            g.glp_set_row_name(lp, i_offset, the_metabolite.id)

            lower_bound = upper_bound = the_metabolite._bound
            constraint_sense = sense_dict[the_metabolite._constraint_sense]
            if constraint_sense == 'E':
                bound_type = GLPKConstants.GLP_FX
            elif constraint_sense == 'L':
                bound_type = GLPKConstants.GLP_UP
            elif constraint_sense == 'G':
                bound_type = GLPKConstants.GLP_LO
            elif constraint_sense == 'U':
                bound_type = GLPKConstants.GLP_FR
            elif hasattr(lower_bound, '__iter__'):
                lower_bound, upper_bound = lower_bound[:2]
                bound_type = GLPKConstants.GLP_DB
                

            g.glp_set_row_bnds(lp, i_offset, bound_type,
                               lower_bound, upper_bound)

            [(row_indices.append(i_offset),
              column_indices.append(reaction_to_index[k]),
              constraint_values.append(k._metabolites[the_metabolite]))
             for k in the_metabolite._reaction]

        #Load the constraints into the lp.  Need to use
        #typed arrays.
        number_of_constraints = len(row_indices)
        i_array = g.new_intArray(number_of_constraints)
        j_array = g.new_intArray(number_of_constraints)
        v_array = g.new_doubleArray(number_of_constraints)
        for a, (i, j, v) in enumerate(zip(row_indices,
                                          column_indices,
                                          constraint_values)):
            g.intArray_setitem(i_array, a+1, i)
            g.intArray_setitem(j_array, a+1, j)
            g.doubleArray_setitem(v_array, a+1, v)
        g.glp_load_matrix(lp, number_of_constraints, i_array,
                          j_array, v_array)
        # the following lines often cause memory crashes
        g.delete_intArray(i_array)
        g.delete_intArray(j_array)
        g.delete_doubleArray(v_array)
        

        g.glp_set_obj_name(lp, "z")
        [g.glp_set_obj_coef(lp, k, v)
          for k, v in iteritems(objective_dict)]

        




__solver_class = Problem

def set_parameter(lp, parameter_name, parameter_value):
    lp.set_parameter(parameter_name, parameter_value)


def get_status(lp):
    return lp.get_status()

def format_solution(lp, cobra_model, **kwargs):
    """

    """
    status = get_status(lp)
    if not lp._mip:
        try:
            x = [lp._g.glp_get_col_prim(lp._lp, i + 1)
                 for i in range(len(cobra_model.reactions))]
            x_dict = dict(zip(cobra_model.reactions, x))

            y = [lp._g.glp_get_row_dual(lp._lp, i + 1)
                 for i in range(len(cobra_model.metabolites))]
            y_dict = dict(zip(cobra_model.metabolites, y))
        
            objective_value = lp.objective_value
        except Exception as e:
            print(repr(e))
            y = y_dict = x = x_dict = objective_value = None
            #print status
    else:
        try:
            x = [lp._g.glp_mip_col_val(lp._lp, i + 1)
                 for i in range(len(cobra_model.reactions))]
            x_dict = dict(zip(cobra_model.reactions, x))
            y = y_dict = None
            objective_value = lp.objective_value
        except:
            y = y_dict = x = x_dict = objective_value = None

    return(Solution(objective_value, x=x, x_dict=x_dict, y=y,
                    y_dict=y_dict, status=status))
def create_problem(cobra_model,  **kwargs):
    """Solver-specific method for constructing a solver problem from
    a cobra.Model.  This can be tuned for performance using kwargs


    """
    the_parameters = parameter_defaults
    if kwargs:
        the_parameters = deepcopy(parameter_defaults)
        the_parameters.update(kwargs)
    quadratic_component = the_parameters['quadratic_component']
    new_objective = the_parameters['new_objective']
    if quadratic_component is not None:
        raise Exception('%s cannot solve QPs, try a different solver'%solver_name)
    lp = Problem()        # Create empty problem instance
    lp.create_problem(cobra_model)
    [set_parameter(lp, parameter_mappings[k], v)
     for k, v in iteritems(the_parameters) if k in parameter_mappings]
    return(lp)

def update_problem(lp, cobra_model, **kwargs):
    """
    Assumes that neither Metabolites nor Reaction have been
    added or removed.

    Currently only deals with reaction bounds and objective
    coefficients.

    """
    g = lp._g
    l = lp._lp
    for i, the_reaction in enumerate(cobra_model.reactions):
        lower_bound = float(the_reaction.lower_bound)
        upper_bound = float(the_reaction.upper_bound)
        objective_coefficient = float(the_reaction.objective_coefficient)
        if lower_bound == upper_bound:
            bound_type = GLPKConstants.GLP_FX
        else:
            bound_type = GLPKConstants.GLP_DB
        g.glp_set_col_bnds(l, i + 1, bound_type, lower_bound, upper_bound)
        g.glp_set_obj_coef(l, i + 1, objective_coefficient)

    
def solve_problem(lp, **kwargs):
    """A performance tunable method for updating a model problem file

    """
    #Update parameter settings if provided
    if kwargs:
        [set_parameter(lp, parameter_mappings[k], v)
         for k, v in iteritems(kwargs) if k in parameter_mappings]
    try:
        print_solver_time = kwargs['print_solver_time']
        start_time = time()
    except:
        print_solver_time = False
    lp_method = lp._simplex_parameters.meth
    lp.solve()
    status = get_status(lp)
    if print_solver_time:
        print('optimize time: {:f}'.format(time() - start_time))
    return status

    
def solve(cobra_model, **kwargs):
    """Smart interface to optimization solver functions that will convert
    the cobra_model to a solver object, set the parameters, and try multiple
    methods to get an optimal solution before returning the solver object and
    a cobra.solution (which is attached to cobra_model.solution)

    cobra_model: a cobra.Model

    returns a dict: {'the_problem': solver specific object, 'the_solution':
    cobra.solution for the optimization problem'}
    

    """
    #Start out with default parameters and then modify if
    #new onese are provided
    the_parameters = deepcopy(parameter_defaults)
    if kwargs:
        the_parameters.update(kwargs)
    #Update objectives if they are new.
    error_reporting = the_parameters['error_reporting']
    if 'new_objective' in the_parameters and \
           the_parameters['new_objective'] not in ['update problem', None]:
       from ..flux_analysis.objective import update_objective
       update_objective(cobra_model, the_parameters['new_objective'])
    if 'the_problem' in the_parameters:
        the_problem = the_parameters['the_problem']
    else:
        the_problem = None
    if isinstance(the_problem, __solver_class):
        #Update the problem with the current cobra_model
        lp = the_problem
        update_problem(lp, cobra_model, **the_parameters)
    else:
        #Create a new problem
        lp = create_problem(cobra_model, **the_parameters)
    #Deprecated way for returning a solver problem created from a cobra_model
    #without performing optimization
    if the_problem == 'setup':
        return lp
    ###Try to solve the problem using other methods if the first method doesn't work
    lp_method = the_parameters['lp_method']
    the_methods = [1, 2, 3]
    if lp_method in the_methods:
        the_methods.remove(lp_method)
    #Start with the user specified method
    the_methods.insert(0, lp_method)
    for the_method in the_methods:
        the_parameters['lp_method'] = the_method
        try:
            status = solve_problem(lp, **the_parameters)
        except:
            status = 'failed'
        if status == 'optimal':
            break
    
    the_solution = format_solution(lp, cobra_model)
    if status != 'optimal' and error_reporting:
        print('{:s} failed: {:s}'.format(solver_name, status))
    cobra_model.solution = the_solution
    solution = {'the_problem': lp, 'the_solution': the_solution}
    return solution
