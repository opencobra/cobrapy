#This centralizes some of the common elements that are differently named across solvers.
#These are stored as strings here to prevent problems associated with calling
#solver objects for solver packages that aren't available
from copy import deepcopy
__objective_sense_cplex = {'maximize': 'Cplex.objective.sense.maximize',
                           'minimize': 'Cplex.objective.sense.minimize'}
__objective_sense_cplex_java = {'maximize': 'IloObjectiveSense.Maximize',
                                'minimize': 'IloObjectiveSense.Minimize'}

__objective_sense_glpk = {'maximize': True,
                          'minimize': False}
__objective_sense_glpk_java = {'maximize': 'GLPKConstants.GLP_MAX',
                               'minimize': 'GLPKConstants.GLP_MIN'}
__objective_sense_gurobi = {'maximize': -1,
                            'minimize': 1}
objective_senses = {'cplex': __objective_sense_cplex,
                    'cplex_java': __objective_sense_cplex_java,
                    'glpk': __objective_sense_glpk,
                    'glpk_java': __objective_sense_glpk_java,
                    'gurobi': __objective_sense_gurobi}
default_objective_sense = 'maximize'
#Mappings from solver-specific status values to cobra pie status values
## __status_cplex = "{Cplex.solution.status.MIP_infeasible: 'infeasible', " +\
##                  "Cplex.solution.status.MIP_optimal: 'optimal', " +\
##                  "Cplex.solution.status.MIP_optimal_tolerance: 'optimal'," +\
##                  "Cplex.solution.status.MIP_unbounded:  'unbounded', "+\
##                  "Cplex.solution.status.infeasible: 'infeasible', " +\
##                  "Cplex.solution.status.optimal: 'optimal',  " +\
##                  "Cplex.solution.status.optimal_tolerance: 'optimal', " +\
##                  "Cplex.solution.status.unbounded: 'unbounded', }"
__status_cplex = """{
    'MIP_infeasible': 'infeasible',
    'integer optimal solution': 'optimal',
    'MIP_optimal': 'optimal',
    'MIP_optimal_tolerance': 'optimal',
    'MIP_unbounded':  'unbounded',
    'infeasible': 'infeasible',
    'optimal': 'optimal',
    'optimal_tolerance': 'optimal',
    'unbounded': 'unbounded',
    'integer optimal, tolerance': 'optimal',
    'time limit exceeded': 'time_limit'
}"""

__status_glpk = """{
    'opt': 'optimal',
    'nofeas': 'infeasible',
    'unbnd': 'unbounded'
}"""
__status_glpk_java = "{GLPKConstants.GLP_OPT: 'optimal', GLPKConstants.GLP_FEAS: 'feasible', GLPKConstants.GLP_INFEAS: 'infeasible', GLPKConstants.GLP_NOFEAS: 'infeasible', GLPKConstants.GLP_UNBND: 'unbounded', GLPKConstants.GLP_UNDEF: 'undefined'}"
__status_gurobi = "{GRB.OPTIMAL: 'optimal', GRB.INFEASIBLE: 'infeasible', GRB.UNBOUNDED: 'unbounded', GRB.TIME_LIMIT: 'time_limit'}"

status_dict = {'cplex': __status_cplex,
               'glpk': __status_glpk,
               'glpk_java': __status_glpk_java,
               'gurobi': __status_gurobi}

#Mappings from solver-specific variable kinds to cobra pie
__kind_cplex = "{'continuous': Cplex.variables.type.continuous, 'integer': Cplex.variables.type.integer}"
__kind_cplex_java = "{'continuous':  IloNumVarType.Float, 'integer': IloNumVarType.Int}"
__kind_glpk = "{'continuous': float, 'integer': int}"
__kind_glpk_java = "{'binary': GLPKConstants.GLP_BV, 'continuous': GLPKConstants.GLP_CV, 'integer': GLPKConstants.GLP_IV}"
__kind_gurobi = "{'continuous': GRB.CONTINUOUS, 'integer': GRB.INTEGER}"

variable_kind_dict = {'cplex': __kind_cplex,
                      'cplex_java': __kind_cplex_java,
                      'glpk': __kind_glpk,
                      'glpk_java': __kind_glpk_java, 
                      'gurobi': __kind_gurobi}

#Mappings from solver-specific constraint senses to cobra pie
sense_dict = {'cplex': "{'E': 'E', 'L': 'L', 'G': 'G'}",
              'glpk': "{'E': 'E', 'L': 'L', 'G': 'G'}",
              'gurobi': "{'E': GRB.EQUAL, 'L': GRB.LESS_EQUAL, 'G': GRB.GREATER_EQUAL}"}


#Mappings from cobra pie parameters names to solver specific parameter names
__mappings_cplex = {'lp_method': 'lpmethod',
                    'lp_parallel': 'threads',
                    'threads': 'threads',
                    'objective_sense': 'objective_sense',
                    'time_limit': 'timelimit',
                    'iteration_limit': 'simplex.limits.iterations',
                    'tolerance_barrier': 'barrier.convergetol',
                    'tolerance_feasibility': 'simplex.tolerances.feasibility',
                    'tolerance_markowitz': 'simplex.tolerances.markowitz',
                    'tolerance_optimality': 'simplex.tolerances.optimality',
                    'MIP_gap_abs': 'mip.tolerances.absmipgap',
                    'MIP_gap': 'mip.tolerances.mipgap'}
__mappings_cplex_java = {'lp_method': 'RootAlg',
                         'lp_parallel': 'ParallelMode',
                         'objective_sense': 'objective_sense',
                         'time_limit': 'TiLim',
                         'tolerance_barrier': 'BarEpComp',
                         'tolerance_feasibility': 'EpRHS',
                         'tolerance_markowitz': 'EpMrk',
                         'tolerance_optimality': 'EpOpt'}
__mappings_glpk = {}
__mappings_glpk_java = {'objective_sense': 'objective_sense',
                        'lp_method': 'meth',
                        'output_verbosity': 'msg_lev',
                        'tolerance_dual': 'tol_dj',
                        'tolerance_integer': 'tol_int',
                        'tolerance_optimality': 'tol_bnd'
                        }
__mappings_gurobi = {'log_file': 'LogFile',
                     'lp_method': 'Method',
                     'threads': 'Threads',
                     'objective_sense': 'ModelSense',
                     'output_verbosity': 'OutputFlag',
                     'quadratic_precision': 'Quad',
                     'time_limit': 'TimeLimit',
                     'tolerance_feasibility': 'FeasibilityTol',
                     'tolerance_markowitz': 'MarkowitzTol',
                     'tolerance_optimality': 'OptimalityTol',
                     'iteration_limit': 'IterationLimit',
                     'MIP_gap_abs': 'MIPGapAbs',
                     'MIP_gap': 'MIPGap'}
parameter_mappings = {'cplex': __mappings_cplex,
                      'cplex_java': __mappings_cplex_java,
                      'glpk': __mappings_glpk,
                      'glpk_java': __mappings_glpk_java,
                      'gurobi': __mappings_gurobi}


#Default solver parameters
__common_defaults = {'objective_sense': 'maximize',
                      'tolerance_optimality': 1e-6, 'tolerance_feasibility': 1e-6,
                      'tolerance_integer': 1e-9, 
                      'quadratic_component': None}
                      


__parameters_cplex = deepcopy(__common_defaults)
__parameters_cplex.update({'lp_method': 1,
                           'lp_parallel': 0,
                           'tolerance_barrier': 1.e-8})
__parameters_glpk = deepcopy(__common_defaults)
__parameters_glpk.update({'lp_method': 1})
__parameters_glpk_java = deepcopy(__common_defaults)
__parameters_glpk_java.update({'lp_method': 1,
                               'output_verbosity': 0,
                               'tolerance_dual': 1e-8})
__parameters_gurobi = deepcopy(__common_defaults)
__parameters_gurobi.update({'output_verbosity': 0,
                          'lp_method': 0,
                          'log_file': '',
                          'tolerance_barrier': 1e-8})


parameter_defaults = {'cplex': __parameters_cplex,
                      'glpk': __parameters_glpk,
                      'glpk_java': __parameters_glpk_java,
                      'gurobi': __parameters_gurobi}
