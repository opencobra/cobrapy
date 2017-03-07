# -*- coding: utf-8 -*-

# Solvers are expected to follow the following interface
# create_problem: makes a solver problem object from a cobra.model and
# sets parameters (if possible)

# format_solution: Returns a cobra.Solution object.  This is where one
# should dress the cobra.model with results if desired.

# get_status: converts a solver specific status flag to a cobra pie flag.

# set_parameter: takes solver specific parameter strings and sets them.

# solve: solves the optimization problem.  this is where one should put
# in logic on what to try if the problem
# isn't optimal

# solve_problem: dumb and fast which will set parameters, if provided

# update_problem: changes bounds and linear objective coefficient of the
# solver specific problem file, given the complementary cobra.model

# This attempts to import all working solvers in this directory

from __future__ import absolute_import

import logging
from os import listdir, path
from warnings import warn

LOGGER = logging.getLogger(__name__)

solver_dict = {}
possible_solvers = set()


def add_solver(solver_name, use_name=None):
    """add a solver module to the solvers"""
    exec("from . import " + solver_name)
    solver = eval(solver_name)
    if use_name is None:
        if hasattr(solver, "solver_name"):
            use_name = solver.solver_name
        else:
            use_name = solver_name
    solver_dict[use_name] = eval(solver_name)


for i in listdir(path.dirname(path.abspath(__file__))):
    if i.startswith("_") or i.startswith(".") or i.startswith('legacy'):
        continue
    if i.startswith("parameters"):
        continue
    if i.endswith(".py") or i.endswith(".so") or i.endswith(".pyc") \
            or i.endswith(".pyd"):
        possible_solvers.add(i.split(".")[0])

if "wrappers" in possible_solvers:
    possible_solvers.remove("wrappers")

for solver in possible_solvers:
    LOGGER.debug("adding '%s'...", solver)
    try:
        add_solver(solver)
    except Exception as err:
        LOGGER.debug("addition failed: %s", str(err))
        pass
    else:
        LOGGER.debug("success!")
    del solver

if len(solver_dict) == 0:
    warn("No LP solvers found")

# clean up the namespace
del path, listdir, warn, i, possible_solvers


class SolverNotFound(Exception):
    None


def get_solver_name(mip=False, qp=False):
    """returns a solver name

    raises SolverNotFound if a suitable solver is not found
    """
    if len(solver_dict) == 0:
        raise SolverNotFound("no solvers installed")
    # glpk only does lp, not qp. Gurobi and cplex are better at mip
    mip_order = ["gurobi", "cplex", "mosek", "coin", "cglpk", "glpk"]
    lp_order = ["cglpk", "cplex",  "gurobi", "mosek", "coin", "glpk"]
    qp_order = ["gurobi", "cplex", "mosek"]

    if mip is False and qp is False:
        for solver_name in lp_order:
            if solver_name in solver_dict:
                return solver_name
        # none of them are in the list order - so return the first one
        return list(solver_dict)[0]
    elif qp:  # mip does not yet matter for this determination
        for solver_name in qp_order:
            if solver_name in solver_dict:
                return solver_name
        # see if any solver defines set_quadratic_objective
        for solver_name in solver_dict:
            if hasattr(solver_dict[solver_name], "set_quadratic_objective"):
                return solver_name
        raise SolverNotFound("no qp-capable solver found")
    else:
        for solver_name in mip_order:
            if solver_name in solver_dict:
                return solver_name
        for solver_name in solver_dict:
            if hasattr(solver_dict[solver_name], "_SUPPORTS_MIP"):
                return solver_name
    raise SolverNotFound("no mip-capable solver found")


def optimize(cobra_model, solver=None, **kwargs):
    """Wrapper to optimization solvers

    solver : str
        Name of the LP solver from solver_dict to use. If None is given, the
        default one will be used

    """
    # If the default solver is not installed then use one of the others
    if solver is None:
        qp = "quadratic_component" in kwargs and \
            kwargs["quadratic_component"] is not None
        solver = get_solver_name(qp=qp)

    return solver_dict[solver].solve(cobra_model, **kwargs)
