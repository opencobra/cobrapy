#cobra.solvers.py
#Converts cobra.Model objects into problems for the LP/QP solvers.

#TODO: Speed up problem construction for each of the optimize routines.
from cobra.core.Solution import Solution
from cobra.flux_analysis.objective import update_objective
from time import time
#This centralizes some of the common elements that are differently named across solvers.
#These are stored as strings here to prevent problems associated with calling
#solver objects for solver packages that aren't available

status_dict = {'cplex': "{'MIP_infeasible': 'infeasible', 'MIP_optimal': 'optimal', 'MIP_optimal_tolerance': 'optimal',  'MIP_unbounded': 'unbounded', 'infeasible': 'infeasible', 'optimal': 'optimal', 'optimal_tolerance': 'optimal', 'unbounded': 'unbounded'}",
               'glpk': "{'opt': 'optimal', 'nofeas': 'infeasible', 'unbnd': 'unbounded'}",
               'gurobi': "{GRB.OPTIMAL: 'optimal', GRB.INFEASIBLE: 'infeasible', GRB.UNBOUNDED: 'unbounded'}"}
variable_kind_dict = {'cplex': "{'continuous': Cplex.variables.type.continuous, 'integer': Cplex.variables.type.integer}",
                      'glpk': "{'continuous': float, 'integer': int}",
                      'gurobi': "{'continuous': GRB.CONTINUOUS, 'integer': GRB.INTEGER}"}


            
def _optimize_cplex(cobra_model, new_objective=None, objective_sense='maximize',
                   min_norm=0, the_problem=None, 
                   tolerance_optimality=1e-6, tolerance_feasibility=1e-6, tolerance_integer=1e-9,
                   tolerance_barrier=1e-8,error_reporting=None, 
                   print_solver_time=False, lp_method=1, lp_parallel=0, copy_problem=False,
                   relax_b=None, quadratic_component=None, reuse_basis=True,
                   update_problem_reaction_bounds=True):
    """Uses the ILOG/CPLEX (www.ibm.com/software/integration/optimization/cplex-optimizer/)
    optimizer to perform an optimization on cobra_model for the objective_coefficients in
    cobra_model._objective_coefficients based on the objective sense.

    cobra_model: A cobra.Model object

    new_objective: Reaction, String, or Integer referring to a reaction in
    cobra_model.reactions to set as the objective.  Currently, only supports single
    objective coeffients.  Will expand to include mixed objectives.

    objective_sense: 'maximize' or 'minimize'

    min_norm: not implemented

    the_problem: None or a problem object for the specific solver that can be used to hot
    start the next solution.

    tolerance_optimality: Solver tolerance for optimality.

    tolerance_feasibility: Solver tolerance for feasibility.

    error_reporting: None or True to disable or enable printing errors encountered
    when trying to find the optimal solution.
    
    print_solver_time: False or True.  Indicates if the time to calculate the solution
    should be displayed.

    quadratic_component: None or 
          scipy.sparse.dok of dim(len(cobra_model.reactions),len(cobra_model.reactions))
         If not None:
          Solves quadratic programming problems for cobra_models of the form:
          minimize: 0.5 * x' * quadratic_component * x + cobra_model._objective_coefficients' * x
          such that,
            cobra_model._lower_bounds <= x <= cobra_model._upper_bounds
            cobra_model._S * x (cobra_model._constraint_sense) cobra_model._b
            
    reuse_basis: Boolean.  If True and the_problem is a model object for the solver,
    attempt to hot start the solution.


    update_problem_reaction_bounds: Boolean.  Set to True if you're providing the_problem
    and you've modified reaction bounds on your cobra_model since creating the_problem.  Only
    necessary for CPLEX

    method for linear optimization: 0 = automatic
    1 = primal simplex, 2 = dual simplex, 3 = network simplex,
    4 = barrier, 5 = sifting, 6 = concurrent dual, barrier, and primal
    
    lp.solve() with Salmonella model:
         cold start: 0.05 seconds
         hot start: 0.05 seconds (slow due to copying the LP)

    """
    if relax_b is not None:
        raise Exception('Need to reimplement constraint relaxation')
    from numpy import array, nan, zeros
    from cobra.flux_analysis.objective import update_objective
    from cobra.solvers.legacy import status_dict, variable_kind_dict

    if error_reporting == 'time' or print_solver_time:
        from time import time
        start_time = time()
    try:
        from cplex import Cplex, SparsePair
        variable_kind_dict = eval(variable_kind_dict['cplex'])
        status_dict = eval(status_dict['cplex'])
    except ImportError as e:
        import sys
        if 'wrong architecture' in e[0] and sys.maxsize > 2**32:
            print 'CPLEX python API is not 64-bit.  please contact your IBM representative'
        else:
            print e
    if new_objective and new_objective != 'update problem':
       update_objective(cobra_model, new_objective)
    if the_problem == None or the_problem in ['return', 'setup', 'parallel'] \
           or not isinstance(the_problem, Cplex):
        lp = Cplex()
        #Using the new objects
        #NOTE: This might be slow
        objective_coefficients = []
        lower_bounds = []
        upper_bounds = []
        variable_names = []
        variable_kinds = []
        [(objective_coefficients.append(x.objective_coefficient),
          lower_bounds.append(x.lower_bound),
          upper_bounds.append(x.upper_bound),
          variable_names.append(x.id),
          variable_kinds.append(variable_kind_dict[x.variable_kind]))
         for x in cobra_model.reactions]
        #Cplex decides that the problem is a MIP if variable_kinds are supplied
        #even if there aren't any integers.
        if Cplex.variables.type.integer in variable_kinds:
            lp.variables.add(obj=objective_coefficients,
                             lb=lower_bounds,
                             ub=upper_bounds,
                             names=variable_names,
                             types=variable_kinds)
        else:
            lp.variables.add(obj=objective_coefficients,
                             lb=lower_bounds,
                             ub=upper_bounds,
                             names=variable_names)

        if relax_b:
            raise Exception('need to reimplement relax_b')
            ## range_values = zeros(len(cobra_model.metabolites))
            ## b_values = array([x._bound for x in cobra_model.metabolties])
            ## for the_nonzero in list(b_values.nonzero()[0]):
            ##     range_values[the_nonzero] = -relax_b
        constraint_sense = []
        constraint_names = []
        constraint_limits = []
        [(constraint_sense.append(x._constraint_sense),
          constraint_names.append(x.id),
          constraint_limits.append(x._bound))
         for x in cobra_model.metabolites]
        
        the_linear_expressions = []
        #NOTE: This won't work with metabolites that aren't in any reaction
        for the_metabolite in cobra_model.metabolites:
            variable_list = []
            coefficient_list = []
            for the_reaction in the_metabolite._reaction:
                variable_list.append(the_reaction.id)
                coefficient_list.append(the_reaction._metabolites[the_metabolite])
            the_linear_expressions.append(SparsePair(ind=variable_list,
                                                     val=coefficient_list))
        if quadratic_component is not None:
            if not hasattr(quadratic_component, 'todok'):
                raise Exception('quadratic component must be a scipy.sparse type array')
            quadratic_component_scaled = quadratic_component.todok()

            lp.parameters.emphasis.numerical.set(1)
            for k, v in quadratic_component_scaled.items():
                lp.objective.set_quadratic_coefficients(int(k[0]), int(k[1]), v)
         

        if relax_b:
            lp.linear_constraints.add(lin_expr=the_linear_expressions,
                                      rhs=constraint_limits,
                                      range_values=list(range_values),
                                      senses=constraint_sense,
                                      names=constraint_names)

        else:
            lp.linear_constraints.add(lin_expr=the_linear_expressions,
                                      rhs=constraint_limits,
                                      senses=constraint_sense,
                                      names=constraint_names)

        if error_reporting == 'time':
            print 'setup new problem: ' + repr(time()-start_time)
            start_time = time()
        
        #Set the problem type as cplex doesn't appear to do this correctly
        problem_type = Cplex.problem_type.LP
        if Cplex.variables.type.integer in variable_kinds:
            if quadratic_component is not None:
                problem_type = Cplex.problem_type.MIQP
            else:
                problem_type = Cplex.problem_type.MILP
        elif quadratic_component is not None:
            problem_type = Cplex.problem_type.QP
        lp.set_problem_type(problem_type)

    else:
        if copy_problem:
            lp = Cplex(the_problem)
            if error_reporting == 'time':
                print 'copy problem: ' + repr(time()-start_time)
                start_time = time()

        else:
            lp = the_problem

        if new_objective:
            lp.objective.set_linear([(x.id, float(x.objective_coefficient))
                                     for x in cobra_model.reactions])
            if error_reporting == 'time':
                print 'set lp objective: ' + repr(time()-start_time)
                start_time = time()
        #SPEED THIS UP
        if update_problem_reaction_bounds:
            lp.variables.set_upper_bounds([(x.id, float(x.upper_bound))
                                            for x in cobra_model.reactions])
            lp.variables.set_lower_bounds([(x.id, float(x.lower_bound))
                                            for x in cobra_model.reactions])

        if error_reporting == 'time':
            print 'changed all bounds: ' + repr(time()-start_time)
            start_time = time()

    if objective_sense == 'maximize':
        lp.objective.set_sense(lp.objective.sense.maximize)
    else:
        lp.objective.set_sense(lp.objective.sense.minimize)
    if tolerance_optimality < 1e-10:
        lp.parameters.simplex.perturbation.constant.set(1)
        lp.parameters.simplex.pgradient.set(1)
        lp.parameters.emphasis.memory.set(1)
        #lp.parameters.simplex.tolerances.markowitz.set(.01)
        lp.parameters.advance.set(2)

    lp.parameters.simplex.tolerances.optimality.set(tolerance_optimality)
    lp.parameters.simplex.tolerances.feasibility.set(tolerance_feasibility)


    if lp.get_problem_type() in [Cplex.problem_type.LP,
                                 Cplex.problem_type.MILP]:
        lp.parameters.lpmethod.set(lp_method)
    elif lp.get_problem_type() in [Cplex.problem_type.QP,
                                 Cplex.problem_type.MIQP]:
        lp.parameters.qpmethod.set(lp_method)


    if lp_parallel > 1:
        lp.parameters.threads.set(lp_parallel)
    #lp.parameters.parallel.set(lp_parallel)
    lp.parameters.barrier.convergetol.set(tolerance_barrier)

    if the_problem == 'setup':
        return lp

    if not error_reporting:
        lp.set_results_stream(None)
        lp.set_warning_stream(None)
    if  print_solver_time:
        start_time = time()
    if not isinstance(the_problem, Cplex):
        #TODO: set tolerance
        lp.solve()
        # Solve this LP with the simplex method.  Takes about 0.2 s without hot start
        lp.status = lp.solution.status[lp.solution.get_status()]
        if lp.status in status_dict:
            status = status_dict[lp.status]
        else:
            status = 'failed'
    else:
        if isinstance(the_problem, Cplex) and reuse_basis:
            try:
                the_basis = the_problem.solution.basis.get_basis()
                lp.start.set_basis(the_basis[0],the_basis[1])
                #TODO: Determine whether the primal or dual works best for the
                #problem of interest.  For the ME matrix the primal appears to
                #work best
                lp_method = 1
                lp.parameters.preprocessing.presolve.set(0)
                lp.parameters.lpmethod.set(lp_method) 
            except:
                print 'no basis in the_problem'
        #TODO: set tolerance and time limit
        #lp.parameters.timelimit.set()
        lp.solve()
        #If the solver takes more than 0.1 s with a hot start it is likely stuck
        lp.status = lp.solution.status[lp.solution.get_status()]
        if lp.status in status_dict:
            status = status_dict[lp.status]
        else:
            status = 'failed'
        if status != 'optimal':
            #Cycle through the different solver options, if a solution is not found
            for lp_method in (1, 2, 3, 4, 5, 6):
                lp = optimize_cplex(cobra_model, new_objective=new_objective,
                                    objective_sense=objective_sense,
                                    min_norm=min_norm, the_problem=None,
                                    print_solver_time=print_solver_time,
                                    tolerance_optimality=tolerance_optimality,
                                    tolerance_feasibility=tolerance_feasibility,
                                    lp_method=lp_method,
                                    quadratic_component=quadratic_component)['the_problem']
                lp.status = lp.solution.status[lp.solution.get_status()]
                if lp.status in status_dict:
                    status = status_dict[lp.status]
                else:
                    status = 'failed'
                if status == 'optimal':
                    break
    if error_reporting == 'time':
        print 'solver time: ' + repr(time()-start_time) + ' with method ' + repr(lp_method)
        start_time = time()

    if print_solver_time:
        print 'cplex time: %f'%(time() - start_time)
    #TODO: It might be able to speed this up a little.
    if status == 'optimal':
        objective_value = lp.solution.get_objective_value()
        #This can be sped up a little
        x_dict = dict(zip(lp.variables.get_names(),
                     lp.solution.get_values()))
        x = array(lp.solution.get_values())
        x = x.reshape(x.shape[0],1)
        #MIP's don't have duals
        if lp.get_problem_type() in (Cplex.problem_type.MIQP,
                                     Cplex.problem_type.MILP):

            y = y_dict = None
        else:
            y_dict = dict(zip(lp.linear_constraints.get_names(),
                              lp.solution.get_dual_values()))
            y = array(lp.solution.get_dual_values())
            y = y.reshape(y.shape[0],1)
    else:
        x = y = x_dict = y_dict = objective_value = None
        if error_reporting:
            print 'cplex failed: %s'%lp.status

    cobra_model.solution = the_solution = Solution(objective_value, x=x, x_dict=x_dict,
                            status=status, y=y, y_dict=y_dict)
    solution = {'the_problem': lp, 'the_solution': the_solution}
    return solution    
   



def _optimize_gurobi(cobra_model, new_objective=None, objective_sense='maximize',
                    min_norm=0, the_problem=None,
                    tolerance_optimality=1e-6, tolerance_feasibility=1e-6,
                    tolerance_barrier=None, tolerance_integer=1e-9, error_reporting=None,
                    print_solver_time=False, copy_problem=False, lp_method=0,
                    relax_b=None, quad_precision=False, quadratic_component=None,
                    reuse_basis=True, lp_parallel=None, update_problem_reaction_bounds=True):
    """Uses the gurobi (http://gurobi.com) optimizer to perform an optimization on cobra_model
    for the objective_coefficients in cobra_model._objective_coefficients based
    on objective sense.

    cobra_model: A cobra.Model object

    new_objective: Reaction, String, or Integer referring to a reaction in
    cobra_model.reactions to set as the objective.  Currently, only supports single
    objective coeffients.  Will expand to include mixed objectives.

    objective_sense: 'maximize' or 'minimize'

    min_norm: not implemented

    the_problem: None or a problem object for the specific solver that can be used to hot
    start the next solution.

    tolerance_optimality: Solver tolerance for optimality.

    tolerance_feasibility: Solver tolerance for feasibility.

    quad_precision: Boolean.  Whether or not to used quad precision in calculations

    error_reporting: None or True to disable or enable printing errors encountered
    when trying to find the optimal solution.
    
    print_solver_time: False or True.  Indicates if the time to calculate the solution
    should be displayed.


    quadratic_component: None or 
          scipy.sparse.dok of dim(len(cobra_model.reactions),len(cobra_model.reactions))
         If not None:
          Solves quadratic programming problems for cobra_models of the form:
          minimize: 0.5 * x' * quadratic_component * x + cobra_model._objective_coefficients' * x
          such that,
            cobra_model._lower_bounds <= x <= cobra_model._upper_bounds
            cobra_model._S * x (cobra_model._constraint_sense) cobra_model._b

            NOTE: When solving quadratic problems it may be necessary to disable quad_precision
            and use lp_method = 0 for gurobi.

    reuse_basis: Boolean.  If True and the_problem is a model object for the solver,
    attempt to hot start the solution.

    update_problem_reaction_bounds: Boolean.  Set to True if you're providing the_problem
    and you've modified reaction bounds on your cobra_model since creating the_problem.  Only
    necessary for CPLEX
    
    lp_parallel: Not implemented

    lp.optimize() with Salmonella model:
         cold start: 0.063 seconds
         hot start: 0.057 seconds (Slow due to copying the LP)
         

    """
    if relax_b is not None:
        raise Exception('Need to reimplement constraint relaxation')
    from numpy import array, nan, zeros
    #TODO: speed this up
    if objective_sense == 'maximize':
        objective_sense = -1
    else:
        objective_sense = 1
    from gurobipy import Model, LinExpr, GRB, QuadExpr
    sense_dict = {'E': GRB.EQUAL,
                  'L': GRB.LESS_EQUAL,
                  'G': GRB.GREATER_EQUAL}
    from cobra.flux_analysis.objective import update_objective
    from cobra.solvers.legacy import status_dict, variable_kind_dict

    variable_kind_dict = eval(variable_kind_dict['gurobi'])
    status_dict = eval(status_dict['gurobi'])

    #Update objectives if they are new.
    if new_objective and new_objective != 'update problem':
       update_objective(cobra_model, new_objective)
    #Create a new problem
    if not the_problem or the_problem in ['return', 'setup'] or \
           not isinstance(the_problem, Model):
        lp = Model("cobra")
        lp.Params.OutputFlag = 0
        lp.Params.LogFile = ''
        # Create variables
        #TODO:  Speed this up 
        variable_list = [lp.addVar(lb=float(x.lower_bound),
                                   ub=float(x.upper_bound),
                                   obj=objective_sense*float(x.objective_coefficient),
                                   name=x.id,
                                   vtype=variable_kind_dict[x.variable_kind])
                         for x in cobra_model.reactions]
        reaction_to_variable = dict(zip(cobra_model.reactions,
                                        variable_list))
        # Integrate new variables
        lp.update()
        #Set objective to quadratic program
        if quadratic_component is not None:
            if not hasattr(quadratic_component, 'todok'):
                raise Exception('quadratic component must be a scipy.sparse type array')

            quadratic_objective = QuadExpr()
            for (index_0, index_1), the_value in quadratic_component.todok().items():
                quadratic_objective.addTerms(the_value,
                                       variable_list[index_0],
                                       variable_list[index_1])
            lp.setObjective(quadratic_objective, sense=objective_sense)
        #Constraints are based on mass balance
        #Construct the lin expression lists and then add
        #TODO: Speed this up as it takes about .18 seconds
        #HERE
        for the_metabolite in cobra_model.metabolites:
            constraint_coefficients = []
            constraint_variables = []
            for the_reaction in the_metabolite._reaction:
                constraint_coefficients.append(the_reaction._metabolites[the_metabolite])
                constraint_variables.append(reaction_to_variable[the_reaction])
            #Add the metabolite to the problem
            lp.addConstr(LinExpr(constraint_coefficients, constraint_variables),
                         sense_dict[the_metabolite._constraint_sense.upper()],
                         the_metabolite._bound,
                         the_metabolite.id)
    else:
        #When reusing the basis only assume that the objective coefficients or bounds can change
        if copy_problem:
            lp = the_problem.copy()
        else:
            lp = the_problem
        if not reuse_basis:
            lp.reset()
        for the_variable, the_reaction in zip(lp.getVars(),
                                              cobra_model.reactions):
            the_variable.lb = float(the_reaction.lower_bound)
            the_variable.ub = float(the_reaction.upper_bound)
            the_variable.obj = float(objective_sense*the_reaction.objective_coefficient)

    
    if the_problem == 'setup':
        return lp
    if print_solver_time:
        start_time = time()
    lp.update()
    lp.setParam("FeasibilityTol", tolerance_feasibility)
    lp.setParam("OptimalityTol", tolerance_optimality) 
    if tolerance_barrier:
        lp.setParam("BarConvTol", tolerance_barrier)

    if quad_precision:
            lp.setParam("Quad", 1)
    lp.setParam("Method", lp_method)

    #Different methods to try if lp_method fails
    the_methods = [0, 2, 1]
    if lp_method in the_methods:
        the_methods.remove(lp_method)
    if not isinstance(the_problem, Model):
        lp.optimize()
        if lp.status in status_dict:
            status = status_dict[lp.status]
        else:
            status = 'failed'
        if status != 'optimal':
            #Try to find a solution using a different method
            lp.setParam("MarkowitzTol", 1e-2)
            for lp_method in the_methods:
                lp.setParam("Method", lp_method)
                lp.optimize()
                if status_dict[lp.status] == 'optimal':
                    break
    else:
        lp.setParam("TimeLimit", 0.6)
        lp.optimize()
        lp.setParam("TimeLimit", "default")
        if lp.status in status_dict:
            status = status_dict[lp.status]
        else:
            status = 'failed'
        if status != 'optimal':
            lp.setParam("MarkowitzTol", 1e-2)
            #Try to find a solution using a different method
            for lp_method in the_methods:
                lp.setParam("Method", lp_method)
                lp.optimize()
                if status_dict[lp.status] == 'optimal':
                    break
                
            if status_dict[lp.status] != 'optimal':
                lp = optimize_gurobi(cobra_model, new_objective=new_objective, objective_sense=objective_sense,
                                     min_norm=min_norm, the_problem=None, 
                                     print_solver_time=print_solver_time)['the_problem']


    if print_solver_time:
        print 'optimize time: %f'%(time() - start_time)
    x_dict = {}
    y_dict = {}
    y = None
    if lp.status in status_dict:
        status = status_dict[lp.status]
    else:
        status = 'failed'
    if status == 'optimal':
        objective_value = objective_sense*lp.ObjVal
        [x_dict.update({v.VarName: v.X}) for v in lp.getVars()]
        x = array([x_dict[v.id] for v in cobra_model.reactions])
        if lp.isMIP:
            y = y_dict = None #MIP's don't have duals
        else:
            [y_dict.update({c.ConstrName: c.Pi})
             for c in lp.getConstrs()]
            y = array([y_dict[v.id] for v in cobra_model.metabolites])
    else:
        y = y_dict = x = x_dict = None
        objective_value = None
        if error_reporting:
            print 'gurobi failed: %s'%lp.status  
    cobra_model.solution = the_solution = Solution(objective_value, x=x, x_dict=x_dict,
                            y=y, y_dict=y_dict,
                            status=status)
    solution = {'the_problem': lp, 'the_solution': the_solution}
    return solution


def _optimize_glpk(cobra_model, new_objective=None, objective_sense='maximize',
                  min_norm=0, the_problem=None, 
                  tolerance_optimality=1e-6, tolerance_feasibility=1e-6, tolerance_integer=1e-9,
                  error_reporting=None, print_solver_time=False,
                 lp_method=1, quadratic_component=None,
                  reuse_basis=True,
                  #Not implemented
                  tolerance_barrier=None, lp_parallel=None,
                  copy_problem=None, relax_b=None,update_problem_reaction_bounds=True):
    """Uses the GLPK (www.gnu.org/software/glpk/) optimizer via pyglpk
    (http://www.tfinley.net/software/pyglpk/release.html) to perform an optimization
    on cobra_model for the objective_coefficients in cobra_model._objective_coefficients
    based on the objective sense.

    cobra_model: A cobra.Model object

    new_objective: Reaction, String, or Integer referring to a reaction in
    cobra_model.reactions to set as the objective.  Currently, only supports single
    objective coeffients.  Will expand to include mixed objectives.

    objective_sense: 'maximize' or 'minimize'

    min_norm: not implemented

    the_problem: None or a problem object for the specific solver that can be used to hot
    start the next solution.

    tolerance_optimality: Solver tolerance for optimality.

    tolerance_feasibility: Solver tolerance for feasibility.

    error_reporting: None or True to disable or enable printing errors encountered
    when trying to find the optimal solution.
    
    print_solver_time: False or True.  Indicates if the time to calculate the solution
    should be displayed.

    quadratic_component: None.  GLPK cannot solve quadratic programs at the moment.

    reuse_basis: Boolean.  If True and the_problem is a model object for the solver,
    attempt to hot start the solution.  Currently, only True is available for GLPK

    
    update_problem_reaction_bounds: Boolean.  Set to True if you're providing the_problem
    and you've modified reaction bounds on your cobra_model since creating the_problem.  Only
    necessary for CPLEX
    
    lp.simplex() with Salmonella model:
         cold start: 0.42 seconds
         hot start: 0.0013 seconds
    """

    from numpy import zeros, array, nan
    #TODO: Speed up problem creation
    if hasattr(quadratic_component, 'todok'):
        raise Exception('GLPK cannot solve quadratic programs please '+\
                        'try using the gurobi or cplex solvers')

    from glpk import LPX
    from cobra.flux_analysis.objective import update_objective
    from cobra.solvers.legacy import status_dict, variable_kind_dict
    status_dict = eval(status_dict['glpk'])
    variable_kind_dict = eval(variable_kind_dict['glpk'])

    if new_objective and new_objective != 'update problem':
        update_objective(cobra_model, new_objective)
    #Faster to use these dicts than index lists
    index_to_metabolite = dict(zip(range(len(cobra_model.metabolites)),
                                   cobra_model.metabolites))
    index_to_reaction = dict(zip(range(len(cobra_model.reactions)),
                                 cobra_model.reactions))
    reaction_to_index = dict(zip(index_to_reaction.values(),
                                 index_to_reaction.keys()))
    if the_problem == None or the_problem in ['return', 'setup'] or \
           not isinstance(the_problem, LPX):
        lp = LPX()        # Create empty problem instance
        lp.name = 'cobra'     # Assign symbolic name to problem
        lp.rows.add(len(cobra_model.metabolites))
        lp.cols.add(len(cobra_model.reactions))
        linear_constraints = []
        for r in lp.rows:
            the_metabolite = index_to_metabolite[r.index]
            r.name = the_metabolite.id
            b = float(the_metabolite._bound)
            c = the_metabolite._constraint_sense
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
    else:
        lp = the_problem
        #BUG with changing / unchanging the basis
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
    if objective_sense.lower() == 'maximize':
        lp.obj.maximize = True # Set this as a maximization problem
    else:
        lp.obj.maximize = False
    if the_problem == 'setup':
        return lp
    if  print_solver_time:
        start_time = time()
    the_methods = [1, 2, 3]
    if lp_method in the_methods:
        the_methods.remove(lp_method)
    else:
        lp_method = 1
    if not isinstance(the_problem, LPX):
       if lp.kind == int:
           lp.simplex(tol_bnd=tolerance_optimality, tol_dj=tolerance_optimality, meth=lp_method)  # we first have to solve the LP?
           lp.integer(tol_int=tolerance_integer)
       else:
           lp.simplex(tol_bnd=tolerance_optimality, tol_dj=tolerance_optimality, meth=lp_method)
       # Solve this LP or MIP with the simplex (depending on if integer variables exist).  Takes about 0.35 s without hot start
       if lp.status in status_dict:
           status = status_dict[lp.status]
       else:
           status = 'failed'
       if status != 'optimal':
           for lp_method in the_methods:
               lp.simplex(tol_bnd=tolerance_optimality, tol_dj=tolerance_optimality, meth=lp_method)
               if lp.status == 'opt':
                   if lp.kind == int:
                       lp.integer(tol_int=tolerance_integer)
                   break
    else:
        if lp.kind == int:
            lp.simplex(tol_bnd=tolerance_optimality, tol_dj=tolerance_optimality, meth=lp_method, tm_lim=100)  # we first have to solve the LP?
            lp.integer(tol_int=tolerance_integer)
        else:
            lp.simplex(tol_bnd=tolerance_optimality, tol_dj=tolerance_optimality, meth=lp_method, tm_lim=100)
       
        #If the solver takes more than 0.1 s with a hot start it is likely stuck
        if lp.status in status_dict:
            status = status_dict[lp.status]
        else:
            status = 'failed'
        if status != 'optimal':
            if lp.kind == int:
               lp.simplex(tol_bnd=tolerance_optimality, tol_dj=tolerance_optimality, meth=lp_method)  # we first have to solve the LP?
               lp.integer(tol_int=tolerance_integer)
            else:
               for lp_method in the_methods:
                   lp.simplex(tol_bnd=tolerance_optimality, tol_dj=tolerance_optimality, meth=lp_method)
                   if lp.status == 'opt':
                       if lp.kind == int:
                           lp.integer(tol_int=tolerance_integer) 
                       break
        if lp.status in status_dict:
            status = status_dict[lp.status]
        else:
            status = 'failed'
        if status != 'optimal':
            lp = optimize_glpk(cobra_model, new_objective=new_objective,
                               objective_sense=objective_sense,
                               min_norm=min_norm, the_problem=None,
                               print_solver_time=print_solver_time,
                               tolerance_optimality=tolerance_optimality,
                               tolerance_feasibility=tolerance_feasibility)['the_problem']
            if lp.status == 'opt':
                if lp.kind == int:
                    lp.integer(tol_int=tolerance_integer)

        if lp.status in status_dict:
            status = status_dict[lp.status]
        else:
            status = 'failed'
        if status != 'optimal':
            lp.simplex(tol_bnd=tolerance_optimality, presolve=True, tm_lim=5000)
            if lp.kind == int:
                lp.integer(tol_int=tolerance_integer)

    if print_solver_time:
        print 'simplex time: %f'%(time() - start_time)
    x = []
    y = []
    x_dict = {}
    y_dict = {}
    if lp.status in status_dict:
        status = status_dict[lp.status]
    else:
        status = 'failed'
    if status == 'optimal':
        objective_value = lp.obj.value
        [(x.append(float(c.primal)),
          x_dict.update({c.name:c.primal}))
          for c in lp.cols]
        
        if lp.kind == float:
            #return the duals as well as the primals for LPs
            [(y.append(float(c.dual)),
              y_dict.update({c.name:c.dual}))
             for c in lp.rows]
        else:
            #MIPs don't have duals
            y = y_dict = None
        x = array(x)
            
    else:
        x = y = x_dict = y_dict = objective_value = None
        if error_reporting:
            print 'glpk failed: %s'%lp.status
    cobra_model.solution = the_solution = Solution(objective_value, x=x, x_dict=x_dict,
                            y=y, y_dict=y_dict,
                            status=status)
    solution = {'the_problem': lp, 'the_solution': the_solution}

    return solution



def test_solvers():
    from cPickle import load
    from time import time
    from numpy import round
    from cobra.manipulation import initialize_growth_medium
    from cobra.test import salmonella_pickle
    with open(salmonella_pickle) as in_file:
        cobra_model = load(in_file)
    initialize_growth_medium(cobra_model, 'M9')
    from cobra.manipulation import initialize_growth_medium
    the_growth_rate = 0.48
    solver_dict = {'glpk': _optimize_glpk,
                   'gurobi': _optimize_gurobi,
                   'cplex': _optimize_cplex}
    try:
        import glpk
    except:
        solver_dict.pop('glpk')
    try:
        from gurobipy import Model
    except:
        solver_dict.pop('gurobi')
    try:
        from cplex import Cplex
    except:
        solver_dict.pop('cplex')
    for the_solver, the_function in solver_dict.items():
        print 'testing ' + the_solver
        start_time = time()
        the_result = the_function(cobra_model, the_problem='return', print_solver_time=True)
        the_problem = the_result['the_problem']
        the_solution = the_result['the_solution']
        print '%s cold start: %f'%(the_solver, time() - start_time)
        if round(the_solution.f, 2) != the_growth_rate:
            print 'Simulation failed %f to match expectation %f'%(the_solution.f,
                                                                  the_growth_rate)
        the_solution = the_function(cobra_model, the_problem=the_problem, print_solver_time=True)
        the_problem = the_result['the_problem']
        the_solution = the_result['the_solution']
        print '%s hot start: %f'%(the_solver, time() - start_time)
        if round(the_solution.f, 2) != the_growth_rate:
            print 'Simulation failed %f to match expectation %f'%(the_solution.f,
                                                                  the_growth_rate)                                                                 
if __name__ == '__main__':
    from cobra.test import create_test_model
    cobra_model = create_test_model()
    _optimize_cplex(cobra_model, error_reporting='time')
