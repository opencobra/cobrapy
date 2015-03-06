from cylp.cy import CyClpSimplex
from cylp.py.modeling.CyLPModel import CyLPArray
from cylp.cy.CyCoinPackedMatrix import CyCoinPackedMatrix

solver_name = "coin"
_status_translation = {"primal infeasible": "infeasible"}

_SUPPORTS_MILP = False


def create_problem(cobra_model, objective_sense="maximize", **kwargs):
    m = cobra_model.to_array_based_model()
    lp = CyClpSimplex()
    v = lp.addVariable("v", len(m.reactions))
    for i, rxn in enumerate(m.reactions):
        if rxn.variable_kind == "integer":
            lp.setInteger(v[i])
    S = m.S
    v.lower = CyLPArray(m.lower_bounds)
    v.upper = CyLPArray(m.upper_bounds)
    inf = float("inf")
    cons = zip(m.b, m.constraint_sense)
    b_l = CyLPArray([-inf if s == "L" else b for b, s in cons])
    b_u = CyLPArray([inf if s == "G" else b for b, s in cons])
    lp.addConstraint(b_u >= S * v >= b_l, "b")
    lp.objectiveCoefficients = CyLPArray(m.objective_coefficients)
    set_parameter(lp, "objective_sense", objective_sense)
    set_parameter(lp, "tolerance_feasibility", 1e-9)
    lp.logLevel = 0
    for key, value in kwargs.items():
        set_parameter(lp, key, value)
    return lp


def solve(cobra_model, **kwargs):
    lp = create_problem(cobra_model)
    for key, value in kwargs.items():
        set_parameter(lp, key, value)
    solve_problem(lp)
    return format_solution(lp, cobra_model)


def set_parameter(lp, parameter_name, value):
    if parameter_name == "objective_sense":
        v = str(value).lower()
        if v == "maximize":
            lp.optimizationDirection = "max"
        elif v == "minimize":
            lp.optimizationDirection = "min"
        else:
            raise ValueError("unknown objective sense '%s'" % value)
    elif parameter_name == "tolerance_feasibility":
        lp.primalTolerance = value
    elif parameter_name == "verbose":
        lp.logLevel = value
    else:
        setattr(lp, parameter_name, value)


def solve_problem(lp, **kwargs):
    for key, value in kwargs.items():
        set_parameter(lp, key, value)
    if max(lp.integerInformation):
        return lp.primal()


def format_solution(lp, cobra_model):
    Solution = cobra_model.solution.__class__
    status = get_status(lp)
    if status != "optimal":  # todo handle other possible
        return Solution(None, status=status)
    solution = Solution(lp.objectiveValue, status=status)
    x = lp.primalVariableSolution["v"].tolist()
    solution.x_dict = {r.id: x[i] for i, r in enumerate(cobra_model.reactions)}
    solution.x = x
    # TODO handle y

    return solution


def get_status(lp):
    status = lp.getStatusString()
    return _status_translation.get(status, status)


def get_objective_value(lp):
    return lp.objectiveValue


def change_variable_bounds(lp, index, lower_bound, upper_bound):
    lp.variablesLower[index] = lower_bound
    lp.variablesUpper[index] = upper_bound


def change_coefficient(lp, met_index, rxn_index, value):
    S = lp.coefMatrix
    S[met_index, rxn_index] = value
    lp.coefMatrix = S


def change_variable_objective(lp, index, value):
    lp.setObjectiveCoefficient(index, value)
