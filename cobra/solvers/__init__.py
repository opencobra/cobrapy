# attempt to import all working solvers in this directory
from os import name as __name
from sys import modules as __modules
from warnings import warn

__legacy_solver = False 
solver_dict = {}
if __legacy_solver:
    from .legacy import _optimize_glpk, _optimize_gurobi, _optimize_cplex
    package_dict = {'glpk': 'from glpk import LPX',
                    'cplex': 'from cplex import Cplex',
                    'gurobi': 'from gurobipy import Model'}
    if __name == 'java':
        from .legacy_jython import _optimize_glpk
        package_dict['glpk'] = 'from org.gnu.glpk import GLPK'
        package_dict['gurobi'] = 'from gurobi import GRBModel'


    solver_dict = {'glpk': _optimize_glpk,
                   'gurobi': _optimize_gurobi,
                   'cplex': _optimize_cplex}

    for solver_name, solver_import in package_dict.iteritems():
        try:
            exec(solver_import)
        except Exception, e:
            #print e
            solver_dict.pop(solver_name)
else:
    from os import listdir as _listdir
    from os import path as _path
    for i in _listdir(_path.split(_path.abspath(__file__))[0]):
        if i.startswith("_") or i.startswith(".") or i == 'legacy.py':
            continue
        if not i.endswith(".py"):
            continue
        try:
            m = i.strip(".py")
            exec("from . import %s" % m)
            solver_name = m
            if solver_name.endswith('_solver'):
                solver_name = solver_name[:-len('_solver')]
            solver_dict[solver_name] = eval(m)
        except Exception, e:
            pass
    del _path
    del _listdir
    del i
    m = None
    del m

def optimize(cobra_model, solver='glpk', error_reporting=True, **kwargs):
    """Wrapper to optimization solvers


    """
    #If the default solver is not installed then use one of the others
    try:
        solver_function = solver_dict[solver]
    except:
        try:
            solver, solver_function = solver_dict.items()[0]
        except:
            raise Exception("It appears that you do not have one of the supported solvers "+\
                            "(glpk, gurobi, or cplex) installed")
    the_solution = None
    if __legacy_solver:
        def solve_problem(solver_function, kwargs):
            return solver_function(cobra_model, **kwargs)
        try:
            the_solution = solve_problem(solver_function, kwargs)
        except Exception, e:
            if error_reporting:
                print e
                print '%s did not work'%solver
            solver_keys = solver_dict.keys()
            solver_keys.remove(solver)
            for solver in solver_keys:
                solver_function = solver_dict[solver]
                try:
                    if error_reporting:
                        print "now trying %s"%solver
                    the_solution = solve_problem(solver_function, kwargs)
                    break
                except Exception, e:
                    if error_reporting:
                        print e
                        print '%s did not work'%solver
                    continue

    else:
        the_solution = solver_function.solve(cobra_model, **kwargs)
        #raise Exception("New style solvers not yet fully implemented")


    #Add the solution to the model.
    #if the_solution is None:
    #   return(the_solution)
    #else:
    return(the_solution['the_problem'])

del __name
