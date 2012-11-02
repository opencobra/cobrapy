##cobra.solvers.glpk_solver
#This script provides wrappers for libglpk-java 1.0.22 and pyglpk 0.3
from warnings import warn
from os import name as __name
from copy import deepcopy
###solver specific parameters
from .parameters import status_dict, variable_kind_dict, \
     sense_dict, parameter_mappings, parameter_defaults, \
     objective_senses, default_objective_sense

from ..core.Solution import Solution
from ..flux_analysis.objective import update_objective
from time import time
solver_name = 'glpk'
sense_dict = eval(sense_dict[solver_name])
#Functions that are different for java implementation of a solver
if __name == 'java':
    warn("cobra.solvers.glpk_solver isn't mature.  consider using gurobi or cplex")
    from org.gnu.glpk import GLPK, GLPKConstants, glp_smcp, glp_iocp
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
                except Exception, e1:
                    try:
                        setattr(self._mip_parameters, parameter_name,
                                parameter_value)
                    except Exception, e2:
                        if warning:
                            print "Could not set simplex parameter " +\
                                  "%s: %s"%(parameter_name, repr(e1))
                            
                            if self._mip_parameters is not None:
                                print "Could not set mip parameter " +\
                                      "%s: %s"%(parameter_name, repr(e2))
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
              for k, v in objective_dict.iteritems()]

            




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
            except Exception, e:
                print repr(e)
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
         for k, v in the_parameters.iteritems() if k in parameter_mappings]
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
             for k, v in kwargs.iteritems() if k in parameter_mappings]
        try:
            print_solver_time = kwargs['print_solver_time']
            start_time = time()
        except:
            print_solver_time = False
        lp_method = lp._simplex_parameters.meth
        lp.solve()
        status = get_status(lp)
        if print_solver_time:
            print 'optimize time: %f'%(time() - start_time)
        return status
else:
    ##Interface to pyGLPK 0.3 
    from glpk import LPX as GLPK
    __solver_class = GLPK
    objective_senses = objective_senses[solver_name]
    variable_kind_dict = eval(variable_kind_dict[solver_name])
    status_dict = eval(status_dict[solver_name])
    parameter_mappings = parameter_mappings[solver_name]
    parameter_defaults = parameter_defaults[solver_name]

    def get_status(lp):
        status = lp.status
        if status in status_dict:
            status = status_dict[status]
        else:
            status = 'failed'
        return status

    def get_objective_value(lp):
        return lp.obj.value

    def format_solution(lp, cobra_model, **kwargs):
        try:
            objective_sign = objective_senses[kwargs['objective_sense']]
        except:
            objective_sign = objective_senses[default_objective_sense]
            
        status = get_status(lp)
        if status == 'optimal':
            objective_value = lp.obj.value
            x = []
            x_dict = {}
            [(x.append(float(c.primal)),
              x_dict.update({c.name:c.primal}))
              for c in lp.cols]

            if lp.kind == float:
                y = []
                y_dict = {}
                #return the duals as well as the primals for LPs
                [(y.append(float(c.dual)),
                  y_dict.update({c.name:c.dual}))
                 for c in lp.rows]
            else:
                #MIPs don't have duals
                y = y_dict = None
            the_solution = Solution(objective_value, x=x, x_dict=x_dict, y=y,
                                    y_dict=y_dict, status=status)
        else:
            the_solution = Solution(None, status=status)
        return(the_solution)

    def set_parameter(lp, parameter_name, parameter_value):
        """with pyglpk the parameters are set during the solve phase, with
        the exception of objective sense.
        
        """
        if parameter_name == 'objective_sense':
            if parameter_value.lower() == 'maximize':
                lp.obj.maximize = True
            elif parameter_value.lower() == 'minimize':
                lp.obj.maximize = False
            else:
                raise ValueError("objective_sense should be 'maximize' or 'minimize'")
        else:
            warn("py glpk solver parameters are set during solve_problem")
        pass

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
        #Faster to use these dicts than index lists
        index_to_metabolite = dict(zip(range(len(cobra_model.metabolites)),
                                       cobra_model.metabolites))
        index_to_reaction = dict(zip(range(len(cobra_model.reactions)),
                                     cobra_model.reactions))
        reaction_to_index = dict(zip(index_to_reaction.values(),
                                     index_to_reaction.keys()))


        lp = __solver_class()        # Create empty problem instance
        lp.name = 'cobra'     # Assign symbolic name to problem
        lp.rows.add(len(cobra_model.metabolites))
        lp.cols.add(len(cobra_model.reactions))
        linear_constraints = []
        for r in lp.rows:
            the_metabolite = index_to_metabolite[r.index]
            r.name = the_metabolite.id
            b = float(the_metabolite._bound)
            c = sense_dict[the_metabolite._constraint_sense]
            if c == 'E':
                r.bounds = b, b     # Set metabolite to steady state levels
            elif c == 'L':
                r.bounds = None, b
            elif c == 'G':
                r.bounds = b, None
            #Add in the linear constraints

            for the_reaction in the_metabolite._reaction:
                reaction_index = reaction_to_index[the_reaction]
                the_coefficient = the_reaction._metabolites[the_metabolite]
                linear_constraints.append((r.index, reaction_index,
                                           the_coefficient))
        #Need to assign lp.matrix after constructing the whole list
        lp.matrix = linear_constraints
        objective_coefficients = []

        for c in lp.cols:
            the_reaction = index_to_reaction[c.index]
            c.name = the_reaction.id           
            the_reaction = index_to_reaction[c.index]
            c.kind = variable_kind_dict[the_reaction.variable_kind]
            c.bounds = the_reaction.lower_bound, the_reaction.upper_bound
            objective_coefficients.append(float(the_reaction.objective_coefficient))
        #Add the new objective coefficients to the problem
        lp.obj[:] = objective_coefficients

        # make sure the objective sense is set in create_problem
        if "objective_sense" in the_parameters:
            set_parameter(lp, "objective_sense", the_parameters["objective_sense"])

        return(lp)


    def update_problem(lp, cobra_model, **kwargs):
        """A performance tunable method for updating a model problem file

        lp: A gurobi problem object

        cobra_model: the cobra.Model corresponding to 'lp'

        """
        #When reusing the basis only assume that the objective coefficients or bounds can change
        #BUG with changing / unchanging the basis
        index_to_metabolite = dict(zip(range(len(cobra_model.metabolites)),
                                       cobra_model.metabolites))
        index_to_reaction = dict(zip(range(len(cobra_model.reactions)),
                                     cobra_model.reactions))
        reaction_to_index = dict(zip(index_to_reaction.values(),
                                     index_to_reaction.keys()))

        try:
            new_objective = kwargs['new_objective']
        except:
            new_objective = None
        if new_objective is not None:
            objective_coefficients = []
            for c in lp.cols:      # Iterate over all rows
                the_reaction = index_to_reaction[c.index]
                c.name = the_reaction.id
                c.bounds = the_reaction.lower_bound, the_reaction.upper_bound
                objective_coefficients.append(float(the_reaction.objective_coefficient))
                c.kind = variable_kind_dict[the_reaction.variable_kind]
            #Add the new objective coefficients to the problem
            lp.obj[:] = objective_coefficients
        else:
            for c in lp.cols:      # Iterate over all rows
                the_reaction = index_to_reaction[c.index]
                c.name = the_reaction.id
                c.bounds = the_reaction.lower_bound, the_reaction.upper_bound
                c.kind = variable_kind_dict[the_reaction.variable_kind]


    ###
    def solve_problem(lp, **kwargs):
        """A performance tunable method for updating a model problem file

        lp: a pyGLPK 0.3 problem

        For pyGLPK it is necessary to provide the following parameters, if they
        are not provided then the default settings will be used: tolerance_optimality,
        tolerance_integer, lp_method, and objective_sense

        """
        if kwargs:
            the_parameters = kwargs
        else:
            the_parameters = {}
        #These need to be provided to solve_problem for pyGLPK because it's not possible
        #AFAIK to set these during problem creation.
        __function_parameters = ['tolerance_optimality', 'tolerance_integer', 'lp_method',
                                 'objective_sense']
        for the_parameter in __function_parameters:
            if the_parameter not in the_parameters:
                the_parameters.update({the_parameter: parameter_defaults[the_parameter]})
               
        try:
            print_solver_time = the_parameters['print_solver_time']
            start_time = time()
        except:
            print_solver_time = False
        tolerance_optimality = the_parameters['tolerance_optimality']
        tolerance_integer = the_parameters['tolerance_integer']
        lp_method = the_parameters['lp_method']
        #[set_parameter(lp, parameter_mappings[k], v)
        # for k, v in kwargs.iteritems() if k in parameter_mappings]

        lp.obj.maximize = objective_senses[the_parameters['objective_sense']] 
        if lp.kind == int:
            #For MILPs, it is faster to solve LP then move to MILP
            lp.simplex(tol_bnd=tolerance_optimality,
                       tol_dj=tolerance_optimality, meth=lp_method)  
            lp.integer(tol_int=tolerance_integer)
        else:
            lp.simplex(tol_bnd=tolerance_optimality, tol_dj=tolerance_optimality,
                       meth=lp_method)
        status = get_status(lp)
        if print_solver_time:
            print 'optimize time: %f'%(time() - start_time)
        return status

    
def solve(cobra_model, **kwargs):
    """Smart interface to optimization solver functions that will convert
    the cobra_model to a solver object, set the parameters, and try multiple
    methods to get an optimal solution before returning the solver object and
    a cobra.Solution (which is attached to cobra_model.solution)

    cobra_model: a cobra.Model

    returns a dict: {'the_problem': solver specific object, 'the_solution':
    cobra.Solution for the optimization problem'}
    

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
        print '%s failed: %s'%(solver_name, status)
    cobra_model.solution = the_solution
    solution = {'the_problem': lp, 'the_solution': the_solution}
    return solution
