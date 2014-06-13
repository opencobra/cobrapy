# attempt to import all working solvers in this directory
from __future__ import absolute_import
from warnings import warn as _warn
from os import name as __name

solver_dict = {}


#TODO: Enforce the solver interface
## create_problem: makes a solver problem object from a cobra.model and
## sets parameters (if possible)

## format_solution: Returns a cobra.Solution object.  This is where one
## should dress the cobra.model with results if desired.

## get_status: converts a solver specific status flag to a cobra pie flag.

## set_parameter: takes solver specific parameter strings and sets them.

## solve: solves the optimization problem.  this is where one should put
## in logic on what to try if the problem
## isn't optimal

## solve_problem: dumb and fast which will set parameters, if provided
##note that for some solvers

## update_problem: changes bounds and linear objective coefficient of the
## solver specific problem file, given the complementary cobra.mod
from os import listdir
from os import path

possible_solvers = set()

def add_solver(solver_name, use_name=None):
    """add a solver module to the solvers"""
    if use_name is None:
        use_name = solver_name
    exec("from . import " + solver_name)
    solver_dict[use_name] = eval(solver_name)


for i in listdir(path.dirname(path.abspath(__file__))):
    if i.startswith("_") or i.startswith(".") or i.startswith('legacy'):
        continue
    if i.startswith("parameters"):
        continue
    if i.endswith(".py") or i.endswith(".so") or i.endswith(".pyc") \
            or i.endswith(".pyd"):
        possible_solvers.add(i.split(".")[0])

for solver in possible_solvers:
    nicer_name = solver[:-7] if solver.endswith("_solver") else solver
    try:
        add_solver(solver, nicer_name)  
    except Exception:
        pass
    del solver, nicer_name

del path, listdir
del i, possible_solvers

if len(solver_dict) == 0:
    _warn("No LP solvers found")

def get_solver_name(mip=False, qp=False):
    """returns a solver name"""
    if len(solver_dict) == 0:
        return None
    # glpk only does lp, not qp. Gurobi and cplex are better at mip
    mip_order = ["gurobi", "cplex", "glpk", "cglpk"]
    lp_order = ["glpk", "cglpk", "gurobi", "cplex"]
    qp_order = ["gurobi", "cplex"]
    qp_incapable = ["cglpk", "glpk"]
    
    if mip is False and qp is False:
        for solver_name in lp_order:
            if solver_name in solver_dict:
                return solver_name
    elif qp:  # mip does not yet matter for this determination
        for solver_name in qp_order:
            if solver_name in solver_dict:
                return solver_name
        for solver_name in solver_dict:
            if solver_name not in qp_incapable:
                _warn("could not verify if %s supports qp" % solver_name)
                return solver_name
        return None  # don't want to return glpk
    else:
        for solver_name in mip_order:
            if solver_name in solver_dict:
                return solver_name
    # return any solver at this point
    return solver_dict.keys()[0]

def optimize(cobra_model, solver=None, **kwargs):
    """Wrapper to optimization solvers

    solver : str
        Name of the LP solver from solver_dict to use. If None is given, the
        default one will be used

    """
    #If the default solver is not installed then use one of the others
    if solver is None:
        solver = get_solver_name()
        if solver is None:
            raise Exception("It appears that you do not have a supported solver")

    solver_function = solver_dict[solver]
    the_solution = solver_function.solve(cobra_model, **kwargs)


    #Add the solution to the model.
    #if the_solution is None:
    #   return(the_solution)
    #else:
    return(the_solution)


del __name
