from six import iteritems

from ..manipulation.modify import convert_to_irreversible, revert_to_reversible
from ..solvers import solver_dict, get_solver_name


def optimize_minimal_flux(cobra_model, already_irreversible=False,
                          fraction_of_optimum=1.0, solver=None,
                          desired_objective_value=None, **optimize_kwargs):
    """Perform basic pFBA (parsimonius FBA) and minimize total flux.

    The function attempts to act as a drop-in replacement for optimize. It
    will make the reaction reversible and perform an optimization, then
    force the objective value to remain the same and minimize the total
    flux. Finally, it will convert the reaction back to the irreversible
    form it was in before. See http://dx.doi.org/10.1038/msb.2010.47

    cobra_model : :class:`~cobra.core.Model` object

    already_irreversible : bool, optional
        By default, the model is converted to an irreversible one.
        However, if the model is already irreversible, this step can be
        skipped

    fraction_of_optimum : float, optional
        Fraction of optimum which must be maintained. The original objective
        reaction is constrained to be greater than maximal_value *
        fraction_of_optimum. By default, this option is specified to be 1.0

    desired_objective_value : float, optional
        A desired objective value for the minimal solution that bypasses the
        initial optimization result.

    solver : string of solver name
        If None is given, the default solver will be used.

    Updates everything in-place, returns model to original state at end.
    """

    if len(cobra_model.objective) > 1:
        raise ValueError('optimize_minimal_flux only supports models with'
                         ' a single objective function')

    if 'objective_sense' in optimize_kwargs:
        if optimize_kwargs['objective_sense'] == 'minimize':
            raise ValueError(
                'Minimization not supported in optimize_minimal_flux')
        optimize_kwargs.pop('objective_sense', None)

    # Convert to irreversible, so all reactions will have a positive flux
    convert_to_irreversible(cobra_model)

    solver = solver_dict[get_solver_name() if solver is None else solver]
    lp = solver.create_problem(cobra_model, **optimize_kwargs)
    if not desired_objective_value:
        solver.solve_problem(lp, objective_sense='maximize')
        status = solver.get_status(lp)
        if status != "optimal":
            revert_to_reversible(cobra_model)
            raise ValueError(
                "pFBA requires optimal solution status, not {}".format(status))
        desired_objective_value = solver.get_objective_value(lp)

    for i, reaction in enumerate(cobra_model.reactions):

        if reaction.objective_coefficient != 0:
            # Enforce a certain fraction of the original objective
            target = (desired_objective_value * fraction_of_optimum /
                      reaction.objective_coefficient)
            solver.change_variable_bounds(lp, i, target, reaction.upper_bound)

        # Minimize all reaction fluxes (including objective?)
        solver.change_variable_objective(lp, i, 1)

    solver.solve_problem(lp, objective_sense='minimize', **optimize_kwargs)
    solution = solver.format_solution(lp, cobra_model)

    # Return the model to its original state
    cobra_model.solution = solution
    revert_to_reversible(cobra_model)

    if solution.status == "optimal":
        cobra_model.solution.f = sum([coeff * reaction.x for reaction, coeff in
                                      iteritems(cobra_model.objective)])

    return solution
