"""bindings to the glpk solver through pyglpk"""
# because both this module and the library are named glpk
from __future__ import absolute_import
from glpk import LPX

from ..core.Solution import Solution


# mappers from cobra representation to glpk
variable_kind_dict = {
    'continuous': float,
    'integer': int}
# mappers from glpk representations to cobra
status_dict = {
    'opt': 'optimal',
    'nofeas': 'infeasible',
    'unbnd': 'unbounded'}


def create_problem(cobra_model, objective_sense="maximize", lp=None):
    if lp is None:
        lp = LPX()  # Create empty problem instance
        lp.name = cobra_model.id
        lp.rows.add(len(cobra_model.metabolites))
        lp.cols.add(len(cobra_model.reactions))

    if objective_sense == 'maximize':
        lp.obj.maximize = True
    elif objective_sense == 'minimize':
        lp.obj.maximize = False
    else:
        raise ValueError("objective_sense not 'maximize' or 'minimize'")

    # create metabolites/constraints as rows
    for i, r in enumerate(lp.rows):
        metabolite = cobra_model.metabolites[i]
        r.name = metabolite.id
        b = float(metabolite._bound)
        c = metabolite._constraint_sense
        # constraint sense is set by changing the bounds
        if c == 'E':
            r.bounds = (b, b)
        elif c == 'L':
            r.bounds = (None, b)
        elif c == 'G':
            r.bounds = (b, None)
        else:
            raise ValueError("%s is not a valid constraint_sense" % c)

    # create reactions/variables as columns
    for i, c in enumerate(lp.cols):
        reaction = cobra_model.reactions[i]
        c.name = reaction.id
        c.kind = variable_kind_dict[reaction.variable_kind]
        c.bounds = (reaction.lower_bound, reaction.upper_bound)
        lp.obj[i] = float(reaction.objective_coefficient)

    # create S matrix
    lp.matrix = [(int(i), int(j), c) \
        for (i, j), c in cobra_model.to_array_based_model().S.todok().iteritems()]
    return lp


def solve_problem(lp):
    lp.simplex()
    try:
        status = status_dict[lp.status]
    except:
        status = "unknown error: " + str(lp.status)
    solution = Solution(status)
    solution.status = status
    if status == 'optimal':
        solution.objective_value = lp.obj.value
        solution.x_dict = dict((c.name, c.primal) for c in lp.cols)
        # return the duals as well as the primals for LPs
        if lp.kind == "float":
            solution.y_dict = dict((c.name, c.dual) for c in lp.cols)
        else:
            solution.y_dict = None
    return solution


def solve(cobra_model, objective_sense="maximize", **kwargs):
    return solve_problem(create_problem(cobra_model, objective_sense))
