from six import iteritems
from itertools import chain
import sympy
import logging

from cobra.manipulation.modify import convert_to_irreversible, \
    revert_to_reversible
from cobra.util import linear_reaction_coefficients, set_objective
from cobra.exceptions import SolveError
import cobra.util.solver as sutil

add = sympy.Add._from_args
mul = sympy.Mul._from_args
LOGGER = logging.getLogger(__name__)


def optimize_minimal_flux(model, already_irreversible=False,
                          fraction_of_optimum=1.0, solver=None,
                          desired_objective_value=None, objective=None,
                          **optimize_kwargs):
    """Perform basic pFBA (parsimonious Enzyme Usage Flux Balance Analysi)
    and minimize total flux.

    pFBA [1] adds the minimization of all fluxes the the objective of the
    model. This approach is motivated by the idea that high fluxes have a
    higher enzyme turn-over and that since producing enzymes is costly,
    the cell will try to minimize overall flux while still maximizing the
    original objective function, e.g. the growth rate.

    Parameters
    ----------
    model : cobra.Model
        The model

    already_irreversible : bool, optional
        By default, the model is converted to an irreversible one.
        However, if the model is already irreversible, this step can be
        skipped. Ignored for optlang solvers as not relevant.

    fraction_of_optimum : float, optional
        Fraction of optimum which must be maintained. The original objective
        reaction is constrained to be greater than maximal_value *
        fraction_of_optimum.

    desired_objective_value : float, optional
        A desired objective value for the minimal solution that bypasses the
        initial optimization result. Ignored for optlang solvers, instead,
        define your objective separately and pass using the `objective`
        argument.

    objective : dict or model.solver.interface.Objective
        A desired objective to use during optimization in addition to the
        pFBA objective. Dictionaries (reaction as key, coefficient as value)
        can be used for linear objectives. Not used for non-optlang solvers.

    solver : str, optional
        Name of the solver to be used. If None it will respect the solver set
        in the model (model.solver).

    **optimize_kwargs : additional arguments for legacy solver, optional
        Additional arguments passed to the legacy solver. Ignored for
        optlang solver (those can be configured using
         model.solver.configuration).

    Returns
    -------
    cobra.core.Solution.Solution
        The solution object obtained optimizing the model with the added
        pFBA constraints.

    References
    ----------
    .. [1] Lewis, N. E., Hixson, K. K., Conrad, T. M., Lerman, J. A.,
       Charusanti, P., Polpitiya, A. D., Palsson, B. O. (2010). Omic data
       from evolved E. coli are consistent with computed optimal growth from
       genome-scale models. Molecular Systems Biology, 6,
       390. doi:10.1038/msb.2010.47

    """
    legacy, solver = sutil.choose_solver(model, solver)
    if not legacy:
        return _optimize_minimal_flux_optlang(
            model, objective=objective,
            fraction_of_optimum=fraction_of_optimum)
    else:
        return _optimize_minimal_flux_legacy(
            model, already_irreversible=already_irreversible,
            fraction_of_optimum=fraction_of_optimum,
            solver=solver,
            desired_objective_value=desired_objective_value,
            **optimize_kwargs)


def add_pfba(model, objective=None, fraction_of_optimum=1.0):
    """Add pFBA objective

    Add objective to minimize the summed flux of all reactions to the
    current objective.

    See Also
    -------
    optimize_minimal_flux

    Parameters
    ----------
    model : cobra.core.Model
        The model to add the objective to
    objective :
        An objective to set in combination with the pFBA objective.
    fraction_of_optimum : float
        Fraction of optimum which must be maintained. The original objective
        reaction is constrained to be greater than maximal_value *
        fraction_of_optimum.
    """
    if objective is not None:
        model.objective = objective
    sutil.fix_objective_as_constraint(model, fraction=fraction_of_optimum)
    reaction_variables = ((rxn.forward_variable, rxn.reverse_variable)
                          for rxn in model.reactions)
    variables = chain(*reaction_variables)
    pfba_objective = model.solver.interface.Objective(add(
        [mul((sympy.singleton.S.One, variable))
         for variable in variables]), direction='min', sloppy=True)
    set_objective(model, pfba_objective)


def _optimize_minimal_flux_optlang(model, objective=None,
                                   fraction_of_optimum=1.0):
    """Helper function to perform pFBA with the optlang interface

    Parameters
    ----------
    model : a cobra model
        The model to perform pFBA on
    objective :
        An objective to use in addition to the pFBA constraints.

    Updates everything in-place, returns model to original state at end.
    """
    solution = None
    with model as m:
        add_pfba(m, objective=objective,
                 fraction_of_optimum=fraction_of_optimum)
        try:
            solution = m.optimize(objective_sense='minimize')
        except SolveError as e:
            LOGGER.error("pfba could not determine an optimal solution for "
                         "objective %s" % m.objective)
            raise e
    return solution


def _optimize_minimal_flux_legacy(model, solver, already_irreversible=False,
                                  fraction_of_optimum=1.0,
                                  desired_objective_value=None,
                                  **optimize_kwargs):
    """Perform basic pFBA (parsimonious FBA) and minimize total flux.

    The function attempts to act as a drop-in replacement for optimize. It
    will make the reaction reversible and perform an optimization, then
    force the objective value to remain the same and minimize the total
    flux. Finally, it will convert the reaction back to the irreversible
    form it was in before. See http://dx.doi.org/10.1038/msb.2010.47

    Parameters
    ----------
    model : cobra.Model
        The model

    solver : solver
        The solver object to use

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

    Updates everything in-place, returns model to original state at end.
    """
    objective_reactions = linear_reaction_coefficients(model)
    if len(objective_reactions) > 1:
        raise ValueError('optimize_minimal_flux only supports models with'
                         ' a single objective function')

    if 'objective_sense' in optimize_kwargs:
        if optimize_kwargs['objective_sense'] == 'minimize':
            raise ValueError(
                'Minimization not supported in optimize_minimal_flux')
        optimize_kwargs.pop('objective_sense', None)

    if not already_irreversible:
        convert_to_irreversible(model)

    lp = solver.create_problem(model, **optimize_kwargs)
    if not desired_objective_value:
        solver.solve_problem(lp, objective_sense='maximize')
        status = solver.get_status(lp)
        if status != "optimal":
            revert_to_reversible(model)
            raise ValueError(
                "pFBA requires optimal solution status, not {}".format(status))
        desired_objective_value = solver.get_objective_value(lp)

    for i, reaction in enumerate(model.reactions):

        if reaction.objective_coefficient != 0:
            # Enforce a certain fraction of the original objective
            target = (desired_objective_value * fraction_of_optimum /
                      reaction.objective_coefficient)
            solver.change_variable_bounds(lp, i, target, reaction.upper_bound)

        # Minimize all reaction fluxes (including objective?)
        solver.change_variable_objective(lp, i, 1)

    solver.solve_problem(lp, objective_sense='minimize', **optimize_kwargs)
    solution = solver.format_solution(lp, model)

    # Return the model to its original state
    model.solution = solution
    revert_to_reversible(model)

    if solution.status == "optimal":
        model.solution.f = sum([coeff * reaction.x for reaction, coeff in
                                iteritems(objective_reactions)])

    return solution
