"""
Wrappers for solvers with an object oriented interface. This creates
functions to call the objects' functions.

The create_problem and solve functions are not included because they
are classmethods. They should be included by specifying
create_problem = PROBLEM_CLASS.create_problem
where PROBLEM_CLASS is the solver class (i.e. GLP, esolver, etc.)
"""


def set_objective_sense(lp, objective_sense="maximize"):
    return lp.set_objective_sense(lp, objective_sense=objective_sense)


def change_variable_bounds(lp, *args):
    return lp.change_variable_bounds(*args)


def change_variable_objective(lp, *args):
    return lp.change_variable_objective(*args)


def change_coefficient(lp, *args):
    return lp.change_coefficient(*args)


def set_parameter(lp, parameter_name, value):
    return lp.change_parameter(parameter_name, value)


def solve_problem(lp, **kwargs):
    return lp.solve_problem(**kwargs)


def get_status(lp):
    return lp.get_status()


def get_objective_value(lp):
    return lp.get_objective_value()


def format_solution(lp, cobra_model):
    return lp.format_solution(cobra_model)
