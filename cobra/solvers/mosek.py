# -*- coding: utf-8 -*-
from __future__ import absolute_import, print_function

import mosek
from six import iteritems, string_types
from six.moves import zip

from cobra.core.solution import LegacySolution

env = mosek.Env()

# make sure the mosek environment works
test = env.Task(0, 0)
test.optimize()
del test


solver_name = "mosek"
__mosek_version__ = ".".join(str(i) for i in mosek.getversion())
_SUPPORTS_MILP = True

status_dict = {
    mosek.solsta.dual_infeas_cer: 'infeasible',
    mosek.solsta.prim_infeas_cer: 'infeasible',
    mosek.solsta.optimal: 'optimal',
    mosek.solsta.integer_optimal: 'optimal'}

param_enums = {
    "bi_clean_optimizer": mosek.optimizertype,
    "cache_license": mosek.onoffkey,
    "check_convexity": mosek.checkconvexitytype,
    "compress_statfile": mosek.onoffkey,
    "feasrepair_optimize": mosek.feasrepairtype,
    "infeas_generic_names": mosek.onoffkey,
    "infeas_prefer_primal": mosek.onoffkey,
    "infeas_report_auto": mosek.onoffkey,
    "license_allow_overuse": mosek.onoffkey,
    "license_wait": mosek.onoffkey,
    "mio_branch_dir": mosek.branchdir,
    "mio_branch_priorities_use": mosek.onoffkey,
    "mio_construct_sol": mosek.onoffkey,
    "mio_cont_sol": mosek.miocontsoltype,
    "mio_use_multithreaded_optimizer": mosek.onoffkey,
    "presolve_lindep_use": mosek.onoffkey,
    "presolve_use": mosek.presolvemode,
    "sim_basis_factor_use": mosek.onoffkey,
    "sim_degen": mosek.simdegen
}

param_aliases = {
    "tolerance_feasibility": "basis_tol_x",
    "time_limit": "optimizer_max_time",
    "mip_gap_abs": "mio_tol_abs_gap",
    "mip_gap": "mio_tol_rel_gap",
    "tolerance_integer": "mio_tol_abs_relax_int"
    }


def _verbose_printer(text):
    print(text)


def create_problem(cobra_model, objective_sense="maximize",
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
    lp.putdouparam(mosek.dparam.intpnt_tol_rel_gap, 1e-6)
    lp.putdouparam(mosek.dparam.presolve_tol_aij, 1e-15)
    lp.putdouparam(mosek.dparam.presolve_tol_abs_lindep, 0.)
    lp.putdouparam(mosek.dparam.presolve_tol_s, 0.)
    lp.putdouparam(mosek.dparam.presolve_tol_x, 1e-10)
    lp.putintparam(mosek.iparam.concurrent_priority_intpnt, 0)
    lp.putintparam(mosek.iparam.concurrent_num_optimizers, 1)
    # add reactions/variables
    lp.appendvars(len(cobra_model.reactions))
    lp.putvarboundlist(
        rxn_indexes,
        (mosek.boundkey.ra,) * n_rxns,
        [float(i.lower_bound) for i in cobra_model.reactions],
        [float(i.upper_bound) for i in cobra_model.reactions],
    )
    lp.putclist(
        rxn_indexes,
        [float(i.objective_coefficient) for i in cobra_model.reactions])
    integer_variables = [i for i, r in enumerate(cobra_model.reactions)
                         if r.variable_kind == "integer"]
    lp.putvartypelist(
        integer_variables,
        (mosek.variabletype.type_int,) * len(integer_variables))

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
        for metabolite, stoichiometry in iteritems(reaction._metabolites):
            lp.putaij(cobra_model.metabolites.index(metabolite),
                      i, stoichiometry)
    # set user-supplied parameters
    for key, value in iteritems(solver_parameters):
        set_parameter(lp, key, value)
    return lp


def set_objective_sense(lp, objective_sense):
    if objective_sense == "maximize":
        lp.putobjsense(mosek.objsense.maximize)
    elif objective_sense == "minimize":
        lp.putobjsense(mosek.objsense.minimize)
    else:
        raise ValueError("unknown objective sense '%s'" % objective_sense)


def set_parameter(lp, parameter_name, parameter_value):
    parameter_name = parameter_name.lower()
    parameter_name = param_aliases.get(parameter_name, parameter_name)
    if parameter_name == "verbose":
        if parameter_value:
            for streamtype in mosek.streamtype.values:
                lp.set_Stream(streamtype, _verbose_printer)
        else:
            for streamtype in mosek.streamtype.values:
                lp.set_Stream(streamtype, lambda x: None)
    elif parameter_name == "objective_sense":
        set_objective_sense(lp, parameter_value)
    elif parameter_name == "quadratic_component":
        if parameter_value is not None:
            set_quadratic_objective(lp, parameter_value)
    elif hasattr(mosek.dparam, parameter_name):
        lp.putdouparam(getattr(mosek.dparam, parameter_name), parameter_value)
    # Integer parameter section
    # Many of these are enumerations. Therefore the value should be converted
    # to the appropriate enum one.
    elif isinstance(parameter_value, string_types) and \
            parameter_name in param_enums:
        lp.putintparam(getattr(mosek.iparam, parameter_name),
                       getattr(param_enums[parameter_name], parameter_value))
    elif isinstance(parameter_value, bool) and parameter_name in param_enums:
        if parameter_value:
            set_parameter(lp, parameter_name, "on")
        else:
            set_parameter(lp, parameter_name, "off")
    elif hasattr(mosek.iparam, parameter_name):
        lp.putintparam(getattr(mosek.iparam, parameter_name), parameter_value)
    else:
        raise ValueError("unknown parameter '%s'" % parameter_name)


def solve_problem(lp, **solver_parameters):
    for key, value in iteritems(solver_parameters):
        set_parameter(lp, key, value)
    lp.optimize()
    return get_status(lp)


def solve(cobra_model, **solver_parameters):
    lp = create_problem(cobra_model, **solver_parameters)
    solve_problem(lp)
    return format_solution(lp, cobra_model)


def _get_soltype(lp):
    """get the solution type

    This is bas for LP and itg for MIP and QP

    """
    if lp.getnumintvar() > 0:
        return mosek.soltype.itg
    if lp.getnumqobjnz() > 0:
        return mosek.soltype.itr
    return mosek.soltype.bas


def get_status(lp):
    mosek_status = lp.getsolsta(_get_soltype(lp))
    return status_dict.get(mosek_status, str(mosek_status))


def get_objective_value(lp):
    return lp.getprimalobj(_get_soltype(lp))


def format_solution(lp, cobra_model):
    soltype = _get_soltype(lp)
    mosek_status = lp.getsolsta(soltype)
    status = status_dict.get(mosek_status, str(mosek_status))
    if status != "optimal":
        return LegacySolution(None, status=status)
    solution = LegacySolution(get_objective_value(lp))
    solution.status = status
    x = [0] * len(cobra_model.reactions)
    lp.getxx(soltype, x)
    solution.x = x
    solution.x_dict = {rxn.id: value for rxn, value
                       in zip(cobra_model.reactions, x)}
    if soltype == mosek.soltype.bas:
        y = [0] * len(cobra_model.metabolites)
        lp.gety(mosek.soltype.bas, y)
        solution.y = y
        solution.y_dict = {met.id: value for met, value
                           in zip(cobra_model.metabolites, y)}
    return solution


def change_variable_objective(lp, index, value):
    lp.putcj(index, value)


def change_variable_bounds(lp, index, lower_bound, upper_bound):
    lp.putvarbound(index, mosek.boundkey.ra, lower_bound, upper_bound)


def change_coefficient(lp, met_index, rxn_index, value):
    lp.putaij(met_index, rxn_index, value)


def set_quadratic_objective(lp, quadratic_objective):
    if not hasattr(quadratic_objective, 'todok'):
        raise Exception('quadratic component must be a sparse matrix')
    row_indexes = []
    col_indexes = []
    values = []
    for (index_0, index_1), value in iteritems(quadratic_objective.todok()):
        # specify lower triangular only
        if index_0 >= index_1:
            row_indexes.append(index_0)
            col_indexes.append(index_1)
            values.append(value)
    lp.putqobj(row_indexes, col_indexes, values)
