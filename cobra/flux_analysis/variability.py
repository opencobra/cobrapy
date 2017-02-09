from warnings import warn

from six import iteritems
from cobra.solvers import solver_dict, get_solver_name
import cobra.util.solver as sutil
from sympy import Float


def flux_variability_analysis(cobra_model, reaction_list=None,
                              fraction_of_optimum=1.0, solver=None,
                              **solver_args):
    legacy = False
    if solver is None:
        solver = cobra_model.solver
    elif "optlang-" in solver:
        solver = sutil.interface_to_str(solver)
        solver = sutil.solvers[solver]
    else:
        legacy = True
        solver = solver_dict[solver]

    if reaction_list is None and "the_reactions" in solver_args:
        reaction_list = solver_args.pop("the_reactions")
        warn("the_reactions is deprecated. Please use reaction_list=")
    if reaction_list is None:
        reaction_list = cobra_model.reactions

    if not legacy:
        return _fva_optlang(cobra_model, reaction_list, fraction_of_optimum)
    else:
        return _fva_legacy(cobra_model, reaction_list, fraction_of_optimum,
                           "maximize", solver, **solver_args)


def _fva_legacy(cobra_model, reaction_list, fraction_of_optimum,
                objective_sense, solver, **solver_args):
    """Runs flux variability analysis to find max/min flux values

    cobra_model : :class:`~cobra.core.Model`:

    reaction_list : list of :class:`~cobra.core.Reaction`: or their id's
        The id's for which FVA should be run. If this is None, the bounds
        will be comptued for all reactions in the model.

    fraction_of_optimum : fraction of optimum which must be maintained.
        The original objective reaction is constrained to be greater than
        maximal_value * fraction_of_optimum

    solver : string of solver name
        If None is given, the default solver will be used.

    """
    lp = solver.create_problem(cobra_model)
    solver.solve_problem(lp, objective_sense=objective_sense)
    solution = solver.format_solution(lp, cobra_model)
    if solution.status != "optimal":
        raise ValueError("FVA requires the solution status to be optimal, "
                         "not " + solution.status)
    # set all objective coefficients to 0
    for i, r in enumerate(cobra_model.reactions):
        if r.objective_coefficient != 0:
            f = solution.x_dict[r.id]
            new_bounds = (f * fraction_of_optimum, f)
            solver.change_variable_bounds(lp, i,
                                          min(new_bounds), max(new_bounds))
            solver.change_variable_objective(lp, i, 0.)
    return calculate_lp_variability(lp, solver, cobra_model, reaction_list,
                                    **solver_args)


def calculate_lp_variability(lp, solver, cobra_model, reaction_list,
                             **solver_args):
    """calculate max and min of selected variables in an LP"""
    fva_results = {str(r): {} for r in reaction_list}
    for what in ("minimum", "maximum"):
        sense = "minimize" if what == "minimum" else "maximize"
        for r in reaction_list:
            r_id = str(r)
            i = cobra_model.reactions.index(r_id)
            solver.change_variable_objective(lp, i, 1.)
            solver.solve_problem(lp, objective_sense=sense, **solver_args)
            fva_results[r_id][what] = solver.get_objective_value(lp)
            # revert the problem to how it was before
            solver.change_variable_objective(lp, i, 0.)
    return fva_results


def _fva_optlang(model, reaction_list, fraction):
    """Helper function to perform FVA with the optlang interface.

    Parameters
    ----------
    model : a cobra model
    reaction_list : list of reactions

    Returns
    -------
    dict
        A dictionary containing the results.
    """
    fva_results = {str(r): {} for r in reaction_list}
    prob = model.solver.interface
    with model as m:
        m.solver.optimize()
        if m.solver.status != "optimal":
            raise ValueError("There is no optimal solution "
                             "for the chosen objective!")
        # Add objective as a variable to the model than set to zero
        # This also uses the fraction to create the lower bound for the
        # old objective
        v = prob.Variable("fva_old_objective",
                          lb=fraction * m.solver.objective.value)
        c = prob.Constraint(
            m.solver.objective.expression - v, lb=0.0, ub=0.0,
            name="fva_old_objective_constraint")
        sutil.add_to_solver(m, [v, c])
        model.objective = Float(0.0)  # This will trigger the reset as well
        for what in ("minimum", "maximum"):
            sense = "min" if what == "minimum" else "max"
            for r in reaction_list:
                r_id = str(r)
                r = m.reactions.get_by_id(r_id)
                # The previous objective assignment already triggers a reset
                # so directly update coefs here to not trigger redundant resets
                # in the history manager
                m.solver.objective.set_linear_coefficients(
                    {r.forward_variable: 1, r.reverse_variable: -1})
                m.solver.objective.direction = sense
                m.solver.optimize()
                fva_results[r_id][what] = m.solver.objective.value
                m.solver.objective.set_linear_coefficients(
                    {r.forward_variable: 0, r.reverse_variable: 0})

    return fva_results


def find_blocked_reactions(cobra_model, reaction_list=None,
                           solver=None, zero_cutoff=1e-9,
                           open_exchanges=False, **solver_args):
    """Finds reactions that cannot carry a flux with the current
    exchange reaction settings for cobra_model, using flux variability
    analysis.

    """
    if solver is None:
        solver = get_solver_name()
    if open_exchanges:
        # should not unnecessarily change model
        cobra_model = cobra_model.copy()
        for reaction in cobra_model.reactions:
            if reaction.boundary:
                reaction.lower_bound = min(reaction.lower_bound, -1000)
                reaction.upper_bound = max(reaction.upper_bound, 1000)
    if reaction_list is None:
        reaction_list = cobra_model.reactions
    # limit to reactions which are already 0. If the reactions alread
    # carry flux in this solution, then they can not be blocked.
    solution = solver_dict[solver].solve(cobra_model, **solver_args)
    reaction_list = [i for i in reaction_list
                     if abs(solution.x_dict[i.id]) < zero_cutoff]
    # run fva to find reactions where both max and min are 0
    flux_span_dict = flux_variability_analysis(
        cobra_model, fraction_of_optimum=0., reaction_list=reaction_list,
        solver=solver, **solver_args)
    return [k for k, v in iteritems(flux_span_dict)
            if max(map(abs, v.values())) < zero_cutoff]
