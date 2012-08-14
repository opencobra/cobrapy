from gurobipy import GRB, Model, LinExpr, GRB, QuadExpr, Column

from ..core.Solution import Solution


# mappers from cobra representation to gurobi
sense_dict = {
    'E': GRB.EQUAL,
    'L': GRB.LESS_EQUAL,
    'G': GRB.GREATER_EQUAL}
variable_kind_dict = {
    'continuous': GRB.CONTINUOUS,
    'integer': GRB.INTEGER}
parameter_dict = {
    # "tolerance_feasibility": "FeasibilityTol",
    # "tolerance_optimality": "OptimalityTol",
    "verbose": "OutputFlag",
    "timelimit": "TimeLimit",
    "threads": "Threads"}
default_parameters = {
}

# mappers from gurobi representations to cobra
status_dict = {
    GRB.OPTIMAL: 'optimal',
    GRB.INFEASIBLE: 'infeasible',
    GRB.UNBOUNDED: 'unbounded',
    GRB.TIME_LIMIT: 'time_limit'}


def create_problem(cobra_model, objective_sense="maximize", **parameters):
    lp = Model("cobra")
    lp.Params.OutputFlag = 0

    if objective_sense == 'maximize':
        objective_sign = -1.0
    elif objective_sense == 'minimize':
        objective_sign = 1.0
    else:
        raise ValueError("objective_sense must be 'maximize' or 'minimize'")

    # create metabolites/constraints
    metabolite_constraints = {}
    for metabolite in cobra_model.metabolites:
        metabolite_constraints[metabolite] = \
            lp.addConstr(0.0, sense_dict[metabolite._constraint_sense],
            metabolite._bound, metabolite.id)
    lp.update()

    # create reactions/variables along with S matrix
    for j, reaction in enumerate(cobra_model.reactions):
        constraints = [metabolite_constraints[i] \
            for i in reaction._metabolites]
        stoichiometry = reaction._metabolites.values()
        lp.addVar(
            lb=float(reaction.lower_bound),
            ub=float(reaction.upper_bound),
            obj=objective_sign * reaction.objective_coefficient,
            name=reaction.id,
            vtype=variable_kind_dict[reaction.variable_kind],
            column=Column(stoichiometry, constraints))
    lp.update()
    return lp


def update_problem(cobra_model, lp, objective_sense="maximize"):
    for reaction in cobra_model.reactions:
        pass
    # TODO implement


def solve_problem(lp, **parameters):
    for key, value in parameters.iteritems():
        if key in parameter_dict:
            lp.setParam(parameter_dict[key], value)
        else:
            lp.setParam(key, value)
    lp.optimize()
    try:
        status = status_dict[lp.status]
    except KeyError:
        status = "unkown error: " + str(lp.status)
    solution = Solution(status)
    solution.status = status
    if status in ("optimal", "time_limit"):
        solution.objective_value = lp.ObjVal * -1  # TODO fix sign
        # solution.x = zeros((n,))
        solution.x = []
        solution.x_dict = {}
        for v in lp.getVars():
            solution.x.append(v.X)
            solution.x_dict[v.VarName] = v.X
        if lp.isMIP:
            solution.y_dict = None  # MIP's don't have duals
        else:
            solution.y_dict = \
                dict((c.ConstrName, c.Pi) for c in lp.getConstrs())
    return solution


def solve(cobra_model, objective_sense="maximize", **kwargs):
    gurobi_args = kwargs
    gurobi_object = create_problem(cobra_model, objective_sense, **gurobi_args)
    solution = solve_problem(gurobi_object)
    return solve_problem(create_problem(cobra_model, objective_sense))
