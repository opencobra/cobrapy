from warnings import warn

from ..external.six import iteritems, string_types
from ..core.Metabolite import Metabolite
from ..solvers import solver_dict, get_solver_name


def flux_variability_analysis(cobra_model, reaction_list=None,
                              fraction_of_optimum=1.0, solver=None,
                              objective_sense="maximize", **solver_args):
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
    if reaction_list is None and "the_reactions" in solver_args:
        reaction_list = solver_args.pop("the_reactions")
        from warnings import warn
        warn("the_reactions is deprecated. Please use reaction_list=")
    if reaction_list is None:
        reaction_list = cobra_model.reactions
    else:
        reaction_list = [cobra_model.reactions.get_by_id(i)
                         if isinstance(i, string_types) else i
                         for i in reaction_list]
    solver = solver_dict[get_solver_name() if solver is None else solver]
    lp = solver.create_problem(cobra_model)
    solver.solve_problem(lp, objective_sense=objective_sense)
    solution = solver.format_solution(lp, cobra_model)
    if solution.status != "optimal":
        raise ValueError("FVA requires the solution status to be optimal, not "
                         + solution.status)
    # set all objective coefficients to 0
    for i, r in enumerate(cobra_model.reactions):
        if r.objective_coefficient != 0:
            f = solution.x_dict[r.id]
            new_bounds = (f * fraction_of_optimum, f)
            solver.change_variable_bounds(lp, i,
                                          min(new_bounds), max(new_bounds))
            solver.change_variable_objective(lp, i, 0.)
    # perform fva
    fva_results = {}
    for r in reaction_list:
        i = cobra_model.reactions.index(r)
        fva_results[r.id] = {}
        solver.change_variable_objective(lp, i, 1.)
        solver.solve_problem(lp, objective_sense="maximize", **solver_args)
        fva_results[r.id]["maximum"] = solver.get_objective_value(lp)
        solver.solve_problem(lp, objective_sense="minimize", **solver_args)
        fva_results[r.id]["minimum"] = solver.get_objective_value(lp)
        # revert the problem to how it was before
        solver.change_variable_objective(lp, i, 0.)
    return fva_results


def find_blocked_reactions(cobra_model, reaction_list=None,
                           solver=None, zero_cutoff=1e-9,
                           open_exchanges=False, **kwargs):
    """Finds reactions that cannot carry a flux with the current
    exchange reaction settings for cobra_model, using flux variability
    analysis.

    """
    from warnings import warn
    if solver is None:
        solver = get_solver_name()
    blocked_reactions = []
    if reaction_list is None and "the_reactions" in solver_args:
        reaction_list = solver_args.pop("the_reactions")
        warn("the_reactions is deprecated. Please use reaction_list=")
    if open_exchanges:
        warn('DEPRECATED: Move to using the Reaction.boundary attribute')
        exchange_reactions = [x for x in cobra_model.reactions
                              if x.startswith('EX')]
        for the_reaction in exchange_reactions:
            if the_reaction.lower_bound >= 0:
                the_reaction.lower_bound = -1000
            if the_reaction.upper_bound >= 0:
                the_reaction.upper_bound = 1000
    flux_span_dict = flux_variability_analysis(
        cobra_model, fraction_of_optimum=0., reaction_list=reaction_list,
        solver=solver, **kwargs)
    blocked_reactions = [k for k, v in flux_span_dict.items()
                         if max(map(abs, v.values())) < zero_cutoff]
    return blocked_reactions
