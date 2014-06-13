# distutils: libraries=glpk
# cython: embedsignature=True

from glpk cimport *
from libc.stdlib cimport malloc, free

from tempfile import NamedTemporaryFile as _NamedTemporaryFile  # for pickling
from os import unlink as _unlink

__glpk_version__ = glp_version()

cdef class GLP:
    cdef glp_prob *glp
    cdef glp_smcp parameters
    cdef glp_iocp integer_parameters


    def __cinit__(self):
        self.glp = glp_create_prob()
        glp_set_obj_dir(self.glp, GLP_MAX)  # default is maximize
        glp_init_smcp(&self.parameters)
        glp_init_iocp(&self.integer_parameters)


    def __dealloc__(self):
        glp_delete_prob(self.glp)

    def __init__(self, cobra_model=None):
        cdef int bound_type, index, m, n, n_values, i
        cdef glp_prob *glp
        cdef int *c_rows
        cdef int *c_cols
        cdef double *c_values
        glp_term_out(GLP_OFF)
        
        # initialize parameters
        self.parameters.msg_lev = GLP_MSG_ERR

        if cobra_model is None:
            return
        glp = self.glp
        m = len(cobra_model.metabolites)
        n = len(cobra_model.reactions)
        glp_add_rows(glp, m)
        glp_add_cols(glp, n)

        metabolite_id_to_index = {r.id: i for i, r in enumerate(cobra_model.metabolites, 1)}

        linear_constraint_rows = []
        linear_constraint_cols = []
        linear_constraint_values = []

        # set metabolite/consraint bounds
        for index, metabolite in enumerate(cobra_model.metabolites, 1):
            b = float(metabolite._bound)
            c = metabolite._constraint_sense
            if c == 'E':
                bound_type = GLP_FX  # Set metabolite to steady state levels
            elif c == 'L':
                bound_type = GLP_UP  # x < 2 <==> x has an upper bound of 2
            elif c == 'G':
                bound_type = GLP_LO  # x > 2 <==> x has a lower bound of 2
            else:
                raise Exception("unsupported bound type: %s" % c)
            glp_set_row_bnds(glp, index, bound_type, b, b)
        
        # set reaction/varaiable bounds
        for index, reaction in enumerate(cobra_model.reactions, 1):
            if reaction.variable_kind == "integer":
                if reaction.lower_bound == 0 and reaction.upper_bound == 1:
                    glp_set_col_kind(self.glp, index, GLP_BV)  # binary
                else:
                    glp_set_col_kind(self.glp, index, GLP_IV)
            if reaction.lower_bound == reaction.upper_bound:
                bound_type = GLP_FX
            else:
                bound_type = GLP_DB
            glp_set_col_bnds(glp, index, bound_type,
                             float(reaction.lower_bound), float(reaction.upper_bound))
            glp_set_obj_coef(glp, index, float(reaction.objective_coefficient))

            for metabolite, coefficient in reaction._metabolites.iteritems():
                metabolite_index = metabolite_id_to_index[metabolite.id]
                linear_constraint_rows.append(metabolite_index)
                linear_constraint_cols.append(index)
                linear_constraint_values.append(coefficient)
        
        # set constraint marix
        # first copy the python lists to c arrays
        n_values = len(linear_constraint_rows)
        c_cols = <int *> malloc((n_values + 1) * sizeof(int))
        c_rows = <int *> malloc((n_values + 1) * sizeof(int))
        c_values = <double *> malloc((n_values + 1) * sizeof(double))
        if c_rows is NULL or c_rows is NULL or c_values is NULL:
            raise MemoryError()
        for i in range(n_values):
            c_rows[i + 1] = linear_constraint_rows[i]
            c_cols[i + 1] = linear_constraint_cols[i]
            c_values[i + 1] = float(linear_constraint_values[i])
        # actually set the values
        glp_load_matrix(glp, n_values, c_rows, c_cols, c_values)
        # free the c arrays
        free(c_rows)
        free(c_cols)
        free(c_values)


    #@classmethod  # decorator does not work 
    def create_problem(cls, cobra_model, objective_sense="maximize"):
        problem = cls(cobra_model)
        problem.set_objective_sense(objective_sense)
        return problem
    create_problem = classmethod(create_problem)


    cpdef change_variable_bounds(self, int index, double lower_bound, double upper_bound):
        cdef int bound_type = GLP_DB
        assert index >= 0
        if lower_bound == upper_bound:
            bound_type = GLP_FX
        glp_set_col_bnds(self.glp, index + 1, bound_type, lower_bound, upper_bound)


    def change_coefficient(self, int met_index, int rxn_index, double value):
        cdef int col_length, i
        cdef int *indexes
        cdef double *values
        # glpk uses 1 indexing
        met_index += 1
        rxn_index += 1
        # we first have to get the old column
        col_length = glp_get_mat_col(self.glp, rxn_index, NULL, NULL)
        indexes = <int *> malloc((col_length + 2) * sizeof(int))
        values = <double *> malloc((col_length + 2) * sizeof(double))
        if indexes == NULL or values == NULL:
            raise MemoryError()
        glp_get_mat_col(self.glp, rxn_index, indexes, values)
        # search for duplicate
        for i in range(col_length):
            # if a duplicate exists replace that value and exit
            if indexes[i + 1] == met_index:
                values[i + 1] = value
                glp_set_mat_col(self.glp, rxn_index, col_length, indexes, values)
                return
        # need to add a new entry
        indexes[col_length + 1] = met_index
        values[col_length + 1] = value
        glp_set_mat_col(self.glp, rxn_index, col_length + 1, indexes, values)
        free(indexes)
        free(values)

    def solve_problem(self, **solver_parameters):
        cdef int result
        cdef glp_smcp parameters = self.parameters
        cdef glp_iocp integer_parameters = self.integer_parameters
        cdef glp_prob *glp = self.glp


        for key, value in solver_parameters.items():
            self.set_parameter(key, value)

        # suspend the gil to allow multithreading
        # multithreading must occur with DIFFERENT glp objects
        # calling solve_problem on the same object from 2 different
        # threads at the same time will probably cause problems
        # because glpk itself is not thread safe

        #with nogil:  # we can use this if glpk ever gets thread-safe malloc
        result = glp_simplex(glp, &parameters)
        assert result == 0
        if self.is_mip():
            self.integer_parameters.tm_lim = self.parameters.tm_lim
            self.integer_parameters.msg_lev = self.parameters.msg_lev
            #self.integer_parameters.tol_bnd = self.parameters.tol_bnd
            #self.integer_parameters.tol_piv = self.parameters.tol_piv
            glp_term_out(GLP_OFF)  # prevent verborse MIP output
            #with nogil:
            result = glp_intopt(glp, &integer_parameters)
            glp_term_out(GLP_ON)
            assert result == 0
        return self.get_status()


    def solve(cls, cobra_model, **kwargs):
        problem = cls.create_problem(cobra_model)
        problem.solve_problem(**kwargs)
        solution = problem.format_solution(cobra_model)
        #cobra_model.solution = solution
        #return {"the_problem": problem, "the_solution": solution}
        return solution
    solve = classmethod(solve)


    def get_status(self):
        cdef int result = glp_mip_status(self.glp) if self.is_mip() else glp_get_status(self.glp)
        if result == GLP_OPT:
            return "optimal"
        if result == GLP_FEAS:
            return glp_get_status(self.glp)
        if result == GLP_UNDEF:
            return "undefined"
        if result == GLP_UNBND:
            return "unbounded"
        if result == GLP_NOFEAS:
            return "infeasible"
        return "failed"

    cpdef set_objective_sense(self, objective_sense):
        objective_sense = objective_sense.lower()
        if objective_sense == "maximize":
            glp_set_obj_dir(self.glp, GLP_MAX)
        elif objective_sense == "minimize":
            glp_set_obj_dir(self.glp, GLP_MIN)
        else:
            raise Exception("%s is not a valid objective sense" % objective_sense)

    cpdef set_parameter(self, parameter_name, value):
        """set a solver parameter

        The following parameters are supported
        time_limit: number of seconds
        """
        if parameter_name == "objective_sense":
            self.set_objective_sense(value)
        elif parameter_name == "time_limit":
            self.parameters.tm_lim = 1000 * int(value)
        elif parameter_name == "tolerance_feasibility":
            self.parameters.tol_bnd = float(value)
            self.parameters.tol_dj = float(value)
        elif parameter_name == "tolerance_markowitz":
            self.parameters.tol_piv = float(value)
        elif parameter_name == "tolerance_integer":
            self.integer_parameters.tol_int = float(value)
        elif parameter_name == "mip_gap":
            self.integer_parameters.mpi_gap = float(value)
        elif parameter_name == "output_verbosity":
            if value is False:
                self.parameters.msg_lev = GLP_MSG_ERR
            elif value is None:
                self.parameters.msg_lev = GLP_MSG_OFF
            elif value is True or value == "all":
                self.parameters.msg_lev = GLP_MSG_ALL
            elif value == "normal":
                self.parameters.msg_lev = GLP_MSG_ON


    cpdef get_objective_value(self):
        if self.is_mip():
            return glp_mip_obj_val(self.glp)
        return glp_get_obj_val(self.glp)


    cpdef change_variable_objective(self, int index, double value):
        assert index >= 0
        glp_set_obj_coef(self.glp, index + 1, value)


    cpdef is_mip(self):
        return glp_get_num_int(self.glp) > 0


    def format_solution(self, cobra_model):
        cdef int i, m, n
        cdef glp_prob *glp = self.glp
        Solution = cobra_model.solution.__class__
        status = self.get_status()
        if status != "optimal":  # todo handle other possible
            return Solution(None, status=status)
        solution = Solution(self.get_objective_value(), status=status)
        m = glp_get_num_rows(glp)
        n = glp_get_num_cols(glp)
        x = [0] * n
        if self.is_mip():
            for i in range(1, n + 1):
                    x[i - 1] = glp_mip_col_val(glp, i)
            #x = [glp_mip_col_val(glp, i) for i in range(1, n + 1)]
            solution.x_dict = {rxn.id: x[i] for i, rxn in enumerate(cobra_model.reactions)}
            solution.x = x
        else:
            for i in range(1, n + 1):
                x[i - 1] = glp_get_col_prim(glp, i)
            #x = [glp_get_col_prim(glp, i) for i in range(1, n + 1)]
            solution.x_dict = {rxn.id: x[i] for i, rxn in enumerate(cobra_model.reactions)}
            solution.x = x
            y = [0] * m
            for i in range(1, m + 1):
                y[i - 1] = glp_get_row_dual(glp, i)
            #y = [glp_get_row_dual(glp, i) for i in range(1, m + 1)]
            solution.y_dict = {met.id: y[i] for i, met in enumerate(cobra_model.metabolites)}
            solution.y = y
        return solution


    def __getstate__(self):
        cdef int result
        cdef char *name
        tempfile = _NamedTemporaryFile(mode="r", delete=False)
        name = tempfile.name
        tempfile.close()
        result = glp_write_prob(self.glp, 0, name)
        assert result == 0
        with open(name, "r") as infile:
            state = infile.read()
        _unlink(name)
        return state


    def __reduce__(self):
        return (GLP, (), self.__getstate__())


    def __setstate__(self, state):
        cdef int result
        cdef char *name = NULL
        with _NamedTemporaryFile(mode="w", delete=False) as tempfile:
            name = tempfile.name
            tempfile.write(state)
        result = glp_read_prob(self.glp, 0, name)
        assert result == 0
        _unlink(name)


    def __copy__(self):
        other = GLP()
        glp_copy_prob(other.glp, self.glp, GLP_ON)
        return other

# wrappers for all the functions at the module level
create_problem = GLP.create_problem
def set_objective_sense(lp, objective_sense="maximize"):
    return lp.set_objective_sense(lp, objective_sense=objective_sense)
cpdef change_variable_bounds(lp, int index, double lower_bound, double upper_bound):
    return lp.change_variable_bounds(index, lower_bound, upper_bound)
cpdef change_variable_objective(lp, int index, double value):
    return lp.change_variable_objective(index, value)
cpdef change_coefficient(lp, int met_index, int rxn_index, double value):
    return lp.change_coefficient(met_index, rxn_index, value)
cpdef set_parameter(lp, parameter_name, value):
    return lp.set_parameter(parameter_name, value)
def solve_problem(lp, **kwargs):
    return lp.solve_problem(**kwargs)
cpdef get_status(lp):
    return lp.get_status()
cpdef get_objective_value(lp):
    return lp.get_objective_value()
cpdef format_solution(lp, cobra_model):
    return lp.format_solution(cobra_model)
solve = GLP.solve
