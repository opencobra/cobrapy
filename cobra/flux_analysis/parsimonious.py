# -*- coding: utf-8 -*-

from __future__ import absolute_import

import logging
from warnings import warn
from itertools import chain

import sympy

from cobra.util import solver as sutil
from cobra.manipulation.modify import (
    convert_to_irreversible, revert_to_reversible)
from cobra.util import linear_reaction_coefficients, set_objective
from cobra.core.solution import get_solution

add = sympy.Add._from_args
mul = sympy.Mul._from_args
LOGGER = logging.getLogger(__name__)


def optimize_minimal_flux(*args, **kwargs):
    warn("optimize_minimal_flux has been renamed to pfba", DeprecationWarning)
    return pfba(*args, **kwargs)


def pfba(model, already_irreversible=False,
         fraction_of_optimum=1.0, solver=None,
         desired_objective_value=None, objective=None,
         reactions=None, **optimize_kwargs):
    """Perform basic pFBA (parsimonious Enzyme Usage Flux Balance Analysis)
    to minimize total flux.

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
    solver : str, optional
        Name of the solver to be used. If None it will respect the solver set
        in the model (model.solver).
    desired_objective_value : float, optional
        A desired objective value for the minimal solution that bypasses the
        initial optimization result. Ignored for optlang solvers, instead,
        define your objective separately and pass using the `objective`
        argument.
    objective : dict or model.problem.Objective
        A desired objective to use during optimization in addition to the
        pFBA objective. Dictionaries (reaction as key, coefficient as value)
        can be used for linear objectives. Not used for non-optlang solvers.
    reactions : iterable
        List of reactions or reaction identifiers. Implies `return_frame` to
        be true. Only return fluxes for the given reactions. Faster than
        fetching all fluxes if only a few are needed. Only supported for
        optlang solvers.
    **optimize_kwargs : additional arguments for legacy solver, optional
        Additional arguments passed to the legacy solver. Ignored for
        optlang solver (those can be configured using
         model.solver.configuration).

    Returns
    -------
    cobra.Solution
        The solution object to the optimized model with pFBA constraints added.

    References
    ----------
    .. [1] Lewis, N. E., Hixson, K. K., Conrad, T. M., Lerman, J. A.,
       Charusanti, P., Polpitiya, A. D., Palsson, B. O. (2010). Omic data
       from evolved E. coli are consistent with computed optimal growth from
       genome-scale models. Molecular Systems Biology, 6,
       390. doi:10.1038/msb.2010.47

    """
    legacy, solver = sutil.choose_solver(model, solver)
    if legacy:
        return _pfba_legacy(
            model, already_irreversible=already_irreversible,
            fraction_of_optimum=fraction_of_optimum, solver=solver,
            desired_objective_value=desired_objective_value,
            **optimize_kwargs)
    else:
        model.solver = solver
        return _pfba_optlang(
            model, objective=objective,
            fraction_of_optimum=fraction_of_optimum, reactions=reactions)


def add_pfba(model, objective=None, fraction_of_optimum=1.0):
    """Add pFBA objective

    Add objective to minimize the summed flux of all reactions to the
    current objective.

    See Also
    -------
    pfba

    Parameters
    ----------
    model : cobra.Model
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
    if model.solver.objective.name == '_pfba_objective':
        raise ValueError('model already has pfba objective')
    sutil.fix_objective_as_constraint(model, fraction=fraction_of_optimum)
    reaction_variables = ((rxn.forward_variable, rxn.reverse_variable)
                          for rxn in model.reactions)
    variables = chain(*reaction_variables)
    pfba_objective = model.problem.Objective(add(
        [mul((sympy.singleton.S.One, variable))
         for variable in variables]), direction='min', sloppy=True,
        name="_pfba_objective")
    set_objective(model, pfba_objective)


def _pfba_optlang(model, objective=None, reactions=None,
                  fraction_of_optimum=1.0):
    """Helper function to perform pFBA with the optlang interface

    Not meant to be used directly.

    Parameters
    ----------
    model : a cobra model
        The model to perform pFBA on
    objective :
        An objective to use in addition to the pFBA constraints.
    reactions : iterable
        List of reactions or reaction identifiers.

    Returns
    -------
    cobra.Solution
        The solution to the pFBA optimization.

    Updates everything in-place, returns model to original state at end.
    """
    reactions = model.reactions if reactions is None \
        else model.reactions.get_by_any(reactions)
    with model as m:
        add_pfba(m, objective=objective,
                 fraction_of_optimum=fraction_of_optimum)
        m.solver.optimize()
        solution = get_solution(m, reactions=reactions)
    return solution


def _pfba_legacy(model, solver, already_irreversible=False,
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
        raise ValueError('pfba only supports models with'
                         ' a single objective function')

    if 'objective_sense' in optimize_kwargs:
        if optimize_kwargs['objective_sense'] == 'minimize':
            raise ValueError(
                'Minimization not supported in pfba')
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
    #    model.solution = solution
    revert_to_reversible(model)

    #    if solution.status == "optimal":
    #        model.solution.f = sum([coeff * reaction.x for reaction, coeff in
    #                                iteritems(objective_reactions)])

    return solution
