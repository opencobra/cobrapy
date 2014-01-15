# attempt to import all working solvers in this directory
from __future__ import absolute_import
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
    if i.endswith(".py") or i.endswith(".so"):
        possible_solvers.add(i[:-3])
    if i.endswith(".pyc") or i.endswith(".pyd"):
        possible_solvers.add(i[:-4])

for solver in possible_solvers:
    nicer_name = solver[:-7] if solver.endswith("_solver") else solver
    try:
        add_solver(solver, nicer_name)
    except Exception:
        pass

del path
del listdir
del i
del solver
del nicer_name


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

    the_solution = solver_function.solve(cobra_model, **kwargs)


    #Add the solution to the model.
    #if the_solution is None:
    #   return(the_solution)
    #else:
    return(the_solution['the_problem'])

del __name
