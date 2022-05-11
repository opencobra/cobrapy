"""Provide functions and helpers to obtain minimal growth media."""

import logging
from typing import TYPE_CHECKING, Iterable, Union

import pandas as pd
from optlang.interface import OPTIMAL
from optlang.symbolics import Zero

from .boundary_types import find_boundary_types


if TYPE_CHECKING:
    from cobra import Model, Reaction


logger = logging.getLogger(__name__)


def add_linear_obj(model: "Model") -> None:
    r"""Add a linear version of a minimal medium to the model solver.

    Changes the optimization objective to finding the growth medium
    requiring the smallest total import flux:

    ..math:: minimize \sum_{r_i in import_reactions} |r_i|

    Parameters
    ----------
    model : cobra.Model
        The cobra model to modify.

    """
    coefs = {}
    for rxn in find_boundary_types(model, "exchange"):
        export = len(rxn.reactants) == 1
        if export:
            coefs[rxn.reverse_variable] = 1
        else:
            coefs[rxn.forward_variable] = 1
    model.objective.set_linear_coefficients(coefs)
    model.objective.direction = "min"


def add_mip_obj(model: "Model") -> None:
    """Add a mixed-integer version of a minimal medium to the model.

    Changes the optimization objective to finding the medium with the least
    components:

    minimize size(R) where R part of import_reactions

    Arguments
    ---------
    model : cobra.model
        The model to modify.

    """
    if len(model.variables) > 1e4:
        logger.warning(
            "The MIP version of minimal media is extremely slow for "
            "models that large :("
        )
    exchange_rxns = find_boundary_types(model, "exchange")
    big_m = max(abs(b) for r in exchange_rxns for b in r.bounds)
    prob = model.problem
    coefs = {}
    to_add = []
    for rxn in exchange_rxns:
        export = len(rxn.reactants) == 1
        indicator = prob.Variable("ind_" + rxn.id, lb=0, ub=1, type="binary")
        if export:
            vrv = rxn.reverse_variable
            indicator_const = prob.Constraint(
                vrv - indicator * big_m, ub=0, name="ind_constraint_" + rxn.id
            )
        else:
            vfw = rxn.forward_variable
            indicator_const = prob.Constraint(
                vfw - indicator * big_m, ub=0, name="ind_constraint_" + rxn.id
            )
        to_add.extend([indicator, indicator_const])
        coefs[indicator] = 1
    model.add_cons_vars(to_add)
    model.solver.update()
    model.objective.set_linear_coefficients(coefs)
    model.objective.direction = "min"


def _as_medium(
    exchanges: Iterable["Reaction"], tolerance: float = 1e-6, exports: bool = False
) -> pd.Series:
    """Convert a solution to medium.

    Parameters
    ----------
    exchanges : list of cobra.reaction
        The exchange reactions to consider.
    tolerance : float > 0, optional
        The absolute tolerance for fluxes. Fluxes with an absolute value
        smaller than this number will be ignored (default 1e-6).
    exports : bool, optional
        Whether to return export fluxes as well (default False).

    Returns
    -------
    pandas.Series
        The "medium", meaning all active import fluxes in the solution.

    """
    logger.debug("Formatting medium.")
    medium = pd.Series()
    for rxn in exchanges:
        export = len(rxn.reactants) == 1
        flux = rxn.flux
        if abs(flux) < tolerance:
            continue
        if export:
            medium[rxn.id] = -flux
        elif not export:
            medium[rxn.id] = flux
    if not exports:
        medium = medium[medium > 0]

    return medium


def minimal_medium(
    model: "Model",
    min_objective_value: float = 0.1,
    exports: bool = False,
    minimize_components: Union[bool, int] = False,
    open_exchanges: bool = False,
) -> Union[pd.Series, pd.DataFrame, None]:
    """Find the minimal growth medium for the `model`.

    Finds the minimal growth medium for the `model` which allows for
    model as well as individual growth. Here, a minimal medium can either
    be the medium requiring the smallest total import flux or the medium
    requiring the least components (ergo ingredients), which will be much
    slower due to being a mixed integer problem (MIP).

    Parameters
    ----------
    model : cobra.model
        The model to modify.
    min_objective_value : float > 0 or array-like object, optional
        The minimum growth rate (objective) that has to be achieved
        (default 0.1).
    exports : bool, optional
        Whether to include export fluxes in the returned medium.
        Defaults to False which will only return import fluxes
        (default False).
    minimize_components : bool or int > 0, optional
        Whether to minimize the number of components instead of the total
        import flux. Might be more intuitive if set to True, but may also
        be slow to calculate for large communities. If set to a number `n`
        will return up to `n` alternative solutions all with the same
        number of components (default False).
    open_exchanges : bool or number, optional
        Whether to ignore currently set bounds and make all exchange
        reactions in the `model` possible. If set to a `number`, all
        exchange reactions will be opened with (-`number`, `number`) as
        bounds (default False).

    Returns
    -------
    pandas.Series, pandas.DataFrame or None
        A pandas.Series giving the import flux for each required import
        reaction and (optionally) the associated export fluxes. All exchange
        fluxes are oriented into the import reaction e.g. positive fluxes
        denote imports and negative fluxes exports. If `minimize_components`
        is a number larger than 1, may return a pandas.DataFrame where each
        column is a minimal medium. Returns None, if the minimization is
        infeasible (for instance if min_growth > maximum growth rate).

    Notes
    -----
    Due to numerical issues, the `minimize_components` option will usually
    only minimize the number of "large" import fluxes. Specifically, the
    detection limit is given by ``integrality_tolerance * max_bound`` where
    ``max_bound`` is the largest bound on an import reaction. Thus, if you
    are interested in small import fluxes as well you may have to adjust the
    solver tolerance at first with `model.tolerance = 1e-7` for instance.
    However, this will be *very* slow for large models especially with GLPK.

    """
    exchange_rxns = find_boundary_types(model, "exchange")
    if isinstance(open_exchanges, bool):
        open_bound = 1000
    else:
        open_bound = open_exchanges

    with model as mod:
        if open_exchanges:
            logger.debug(f"Opening exchanges for {len(exchange_rxns)} imports.")
            for rxn in exchange_rxns:
                rxn.bounds = (-open_bound, open_bound)
        logger.debug("Applying objective value constraints.")
        obj_const = mod.problem.Constraint(
            mod.objective.expression,
            lb=min_objective_value,
            name="medium_obj_constraint",
        )
        mod.add_cons_vars([obj_const])
        mod.solver.update()
        mod.objective = Zero
        logger.debug("Adding new media objective.")
        tol = mod.solver.configuration.tolerances.feasibility

        if minimize_components:
            add_mip_obj(mod)
            if isinstance(minimize_components, bool):
                minimize_components = 1
            seen = set()
            best = num_components = mod.slim_optimize()
            if mod.solver.status != OPTIMAL:
                logger.warning("Minimization of medium was infeasible.")
                return None
            exclusion = mod.problem.Constraint(Zero, ub=0)
            mod.add_cons_vars([exclusion])
            mod.solver.update()
            media = []

            for i in range(minimize_components):
                logger.info(f"Finding alternative medium #{(i + 1)}.")
                vars = [mod.variables["ind_" + s] for s in seen]
                if len(seen) > 0:
                    exclusion.set_linear_coefficients(dict.fromkeys(vars, 1))
                    exclusion.ub = best - 1
                num_components = mod.slim_optimize()
                if mod.solver.status != OPTIMAL or num_components > best:
                    if i == 0:
                        logger.warning(
                            "Could not get an optimal solution. "
                            "This is usually due to numerical instability. "
                            "Possible remedies are relaoding the model "
                            "from scratch, switching to a different solver, "
                            "or decreasing the solver tolerance. Please, "
                            "carefully read the note on numerical instability "
                            "in the function documentation."
                        )
                        return None
                    break
                medium = _as_medium(exchange_rxns, tol, exports=exports)
                media.append(medium)
                seen.update(medium[medium > 0].index)
            if len(media) > 1:
                medium = pd.concat(media, axis=1, sort=True).fillna(0.0)
                medium.sort_index(axis=1, inplace=True)
            else:
                medium = media[0]
        else:
            add_linear_obj(mod)
            mod.slim_optimize()
            if mod.solver.status != OPTIMAL:
                logger.warning("Minimization of medium was infeasible.")
                return None
            medium = _as_medium(exchange_rxns, tol, exports=exports)

    return medium
