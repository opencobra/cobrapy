# attempt to import all working solvers in this directory
__legacy_solver = True
solver_dict = {}
if __legacy_solver:
    from .legacy import _optimize_glpk, _optimize_gurobi, _optimize_cplex
    solver_dict = {'glpk': _optimize_glpk,
                   'gurobi': _optimize_gurobi,
                   'cplex': _optimize_cplex}
else:
    from os import listdir as _listdir
    from os import path as _path
    for i in _listdir(_path.split(_path.abspath(__file__))[0]):
        if i.startswith("_") or i.startswith("."):
            continue
        if not i.endswith(".py"):
            continue
        try:
            exec("import .%s" % i.strip(".py"))
            solver_dict[i.strip(".py")] = eval(i.strip(".py"))
        except Exception, e:
            print Exception, e
            pass
    del _path
    del _listdir
    del i
def optimize(cobra_model, solver='glpk', **kwargs):
    """Wrapper to optimization solvers
    

    """
    solver_function = solver_dict[solver]
    the_solution = None
    if __legacy_solver:
        def solve_problem(solver_function, kwargs):
            return solver_function(cobra_model, **kwargs)
        try:
            the_solution = solve_problem(solver_function, kwargs)
        except Exception, e:
            print e
            print '%s did not work'%solver
            solver_keys = solver_dict.keys()
            solver_keys.remove(solver)
            for solver in solver_keys:
                solver_function = solver_dict[solver]
                try:
                    print "now trying %s"%solver
                    the_solution = solve_problem(solver_function, kwargs)
                    break
                except Exception, e:
                    print e
                    print '%s did not work'%solver
                    continue

    else:
        raise Exception("New style solvers not yet fully implemented")


    #Add the solution to the model.
    #Please note that it will be faster to dress the reactions and metabolites
    #with their values within the new style solvers instead of here.
    if the_solution is None:
        cobra_model.solution = the_solution
        return(the_solution)
    else:
        cobra_model.solution = the_solution['the_solution']
        return(the_solution['the_problem'])
