from __future__ import absolute_import

import mosek

env = mosek.Env()

solver_name = "mosek"
__mosek_version__ = ".".join(str(i) for i in mosek.getversion())
_SUPPORTS_MILP = False

status_dict = {
    mosek.solsta.dual_infeas_cer: 'infeasible',
    mosek.solsta.prim_infeas_cer: 'infeasible',
    mosek.solsta.optimal: 'optimal'}


def verbose_printer(text):
        print text


def create_problem(cobra_model, objective_sense="maximize", verbose=False,
                   **solver_parameters):
    n_rxns = len(cobra_model.reactions)
    n_mets = len(cobra_model.metabolites)
    rxn_indexes = range(n_rxns)
    met_indexes = range(n_mets)
    # create lp object and set default parameters
    lp = env.Task(0, 0)
    set_objective_sense(lp, objective_sense)
    lp.putdouparam(mosek.dparam.basis_tol_x, 1e-9)
    lp.putdouparam(mosek.dparam.basis_tol_s, 1e-9)
    lp.putdouparam(mosek.dparam.basis_rel_tol_s, 0.)
    lp.putdouparam(mosek.dparam.simplex_abs_tol_piv, 1e-12)
    lp.putdouparam(mosek.dparam.presolve_tol_aij, 1e-15)
    lp.putdouparam(mosek.dparam.presolve_tol_abs_lindep, 0.)
    lp.putdouparam(mosek.dparam.presolve_tol_s, 0.)
    lp.putdouparam(mosek.dparam.presolve_tol_x, 0.)
    lp.putintparam(mosek.iparam.concurrent_priority_intpnt, 0)
    lp.putintparam(mosek.iparam.concurrent_num_optimizers, 1)
    for key, value in solver_parameters:
        set_parameter(lp, key, value)
    # add reactions/variables
    lp.appendvars(len(cobra_model.reactions))
    lp.putvarboundlist(
        rxn_indexes,
        (mosek.boundkey.ra,) * n_rxns,
        cobra_model.reactions.list_attr("lower_bound"),
        cobra_model.reactions.list_attr("upper_bound"),
    )
    lp.putclist(rxn_indexes,
                cobra_model.reactions.list_attr("objective_coefficient"))

    # add metabolites/constraints
    c_sense_dict = {"E": mosek.boundkey.fx, "L": mosek.boundkey.up,
                    "G": mosek.boundkey.lo}
    c_senses = [c_sense_dict[met._constraint_sense]
                for met in cobra_model.metabolites]
    bounds = cobra_model.metabolites.list_attr("_bound")
    lp.appendcons(len(cobra_model.metabolites))
    lp.putconboundlist(met_indexes, c_senses, bounds, bounds)

    # add in the S matrix
    for i, reaction in enumerate(cobra_model.reactions):
        for metabolite, stoiciometry in reaction._metabolites.iteritems():
            lp.putaij(cobra_model.metabolites.index(metabolite),
                      i, stoiciometry)

    return lp


def set_objective_sense(lp, objective_sense):
    if objective_sense == "maximize":
        lp.putobjsense(mosek.objsense.maximize)
    elif objective_sense == "minimize":
        lp.putobjsense(mosek.objsense.minimize)
    else:
        raise ValueError("unknown objective sense '%s'" % objective_sense)


def set_parameter(lp, parameter_name, parameter_value):
    if parameter_name == "verbose":
        if parameter_value:
            lp.set_Stream(mosek.streamtype.log, verbose_printer)
        else:
            lp.set_Stream(mosek.streamtype.log, lambda x: None)
    elif parameter_name == "objective_sense":
        set_objective_sense(lp, parameter_value)
    elif hasattr(mosek.dparam, parameter_name):
        lp.putdouparam(getattr(mosek.dparam, parameter_name), parameter_value)
    elif hasattr(mosek.iparam, parameter_name):
        lp.putintparam(getattr(mosek.iparam, parameter_name), parameter_value)
    else:
        raise ValueError("unknown parameter '%s'" % parameter_name)


def solve_problem(lp, **solver_parameters):
    for key, value in solver_parameters.items():
        set_parameter(lp, key, value)
    lp.optimize()
    return get_status(lp)


def solve(cobra_model, **solver_parameters):
    lp = create_problem(cobra_model, **solver_parameters)
    solve_problem(lp)
    return format_solution(lp, cobra_model)


def get_status(lp):
    status = lp.getsolsta(mosek.soltype.bas)
    return status_dict.get(status, str(status))


def get_objective_value(lp):
    return lp.getprimalobj(mosek.soltype.bas)


def format_solution(lp, cobra_model):
    status = get_status(lp)
    if status != "optimal":
        return cobra_model.solution.__class__(None, status=status)
    solution = cobra_model.solution.__class__(get_objective_value(lp))
    solution.status = status
    x = [0] * len(cobra_model.reactions)
    y = [0] * len(cobra_model.metabolites)
    lp.getxx(mosek.soltype.bas, x)
    lp.gety(mosek.soltype.bas, y)
    solution.x = x
    solution.y = y
    # TODO: use izip instead of enumerate in dict comprehensions
    solution.x_dict = {rxn.id: x[i] for i, rxn
                       in enumerate(cobra_model.reactions)}
    solution.y_dict = {met.id: y[i] for i, met
                       in enumerate(cobra_model.metabolites)}
    return solution


def change_variable_objective(lp, index, value):
    lp.putcj(index, value)


def change_variable_bounds(lp, index, lower_bound, upper_bound):
    lp.putvarbound(index, mosek.boundkey.ra, lower_bound, upper_bound)


def change_coefficient(lp, met_index, rxn_index, value):
    lp.putaij(met_index, rxn_index, value)
