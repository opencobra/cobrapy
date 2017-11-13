"""Manages functions for growth media analysis and manipulation."""

from collections import Counter
from optlang.symbolics import Zero
from optlang.interface import OPTIMAL
import numpy as np
import pandas as pd
import logging

LOGGER = logging.getLogger(__name__)

default_excludes = ["biosynthesis", "transcription", "replication", "sink",
                    "demand", "DM_", "SN_", "SK_"]
"""A list of sub-strings in reaction IDs that usually indicate that
the reaction is *not* an exchange reaction."""


def external_compartment(model):
    """Find the external compartment in the model.

    Uses a simple heuristic where the external compartment should be the one
    with the most exchange reactions.

    Arguments
    ---------
    model : cobra.Model
        A cobra model.

    Returns
    -------
    str
        The putative external compartment.
    """
    counts = Counter(tuple(r.compartments)[0] for r in model.exchanges)
    return counts.most_common(1)[0][0]


def exchanges(model, ext_compartment=None):
    """Find exchange reactions.

    Arguments
    ---------
    model : cobra.Model
        A cobra model.
    ext_compartment : str or None
        The id for the external compartment. If None will be detected
        automatically.

    Returns
    -------
    list of cobra.reaction
        A list of likely exchange reactions.
    """
    if ext_compartment is None:
        ext_compartment = external_compartment(model)
    exchange_reactions = model.reactions.query(
            lambda r: (r.boundary and not
                       any(ex in r.id for ex in default_excludes) and
                       ext_compartment in r.compartments))
    return exchange_reactions


def add_linear_obj(model):
    """Add a linear version of a minimal medium to the model solver.

    Changes the optimization objective to finding the growth medium requiring
    the smallest total import flux::

        minimize sum |r_i| for r_i in import_reactions

    Arguments
    ---------
    model : cobra.Model
        The model to modify.
    """
    coefs = {}
    for rxn in exchanges(model):
        export = len(rxn.reactants) == 1
        if export:
            coefs[rxn.reverse_variable] = 1
        else:
            coefs[rxn.forward_variable] = 1
    model.objective.set_linear_coefficients(coefs)
    model.objective.direction = "min"


def add_mip_obj(model):
    """Add a mixed-integer version of a minimal medium to the model.

    Changes the optimization objective to finding the medium with the least
    components::

        minimize size(R) where R part of import_reactions

    Arguments
    ---------
    model : cobra.model
        The model to modify.
    """
    if len(model.variables) > 1e4:
        LOGGER.warning("the MIP version of minimal media is extremely slow for"
                       " models that large :(")
    exchange_rxns = exchanges(model)
    M = max(np.max(np.abs(r.bounds)) for r in exchange_rxns)
    prob = model.problem
    coefs = {}
    to_add = []
    for rxn in exchange_rxns:
        export = len(rxn.reactants) == 1
        indicator = prob.Variable("ind_" + rxn.id, lb=0, ub=1, type="binary")
        if export:
            vrv = rxn.reverse_variable
            indicator_const = prob.Constraint(
                vrv - indicator * M, ub=0, name="ind_constraint_" + rxn.id)
        else:
            vfw = rxn.forward_variable
            indicator_const = prob.Constraint(
                vfw - indicator * M, ub=0, name="ind_constraint_" + rxn.id)
        to_add.extend([indicator, indicator_const])
        coefs[indicator] = 1
    model.add_cons_vars(to_add)
    model.solver.update()
    model.objective.set_linear_coefficients(coefs)
    model.objective.direction = "min"


def _as_medium(exchanges, tolerance=1e-6, exports=False):
    """Convert a solution to medium.

    Arguments
    ---------
    exchanges : list of cobra.reaction
        The exchange reactions to consider.
    tolerance : positive double
        The absolute tolerance for fluxes. Fluxes with an absolute value
        smaller than this number will be ignored.
    exports : bool
        Whether to return export fluxes as well.

    Returns
    -------
    pandas.Series
        The "medium", meaning all active import fluxes in the solution.
    """
    LOGGER.info("formatting medium")
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


def minimal_medium(model, min_growth=0.1, exports=False,
                   minimize_components=False, open_exchanges=False):
    """Find the minimal growth medium for the model.

    Finds the minimal growth medium for the model which allows for
    model as well as individual growth. Here, a minimal medium can either
    be the medium requiring the smallest total import flux or the medium
    requiring the least components (ergo ingredients), which will be much
    slower.

    Arguments
    ---------
    model : cobra.model
        The model to modify.
    min_growth : positive float or array-like object.
        The minimum growth rate (objective) that has to be achieved.
    exports : boolean
        Whether to include export fluxes in the returned medium. Defaults to
        False which will only return import fluxes.
    minimize_components : boolean or positive int
        Whether to minimize the number of components instead of the total
        import flux. Might be more intuitive if set to True but may also be
        slow to calculate for large communities. If set to a number `n` will
        return up to `n` alternative solutions all with the same number of
        components.
    open_exchanges : boolean or number
        Whether to ignore currently set bounds and make all exchange reactions
        in the model possible. If set to a number all exchange reactions will
        be opened with (-number, number) as bounds.

    Returns
    -------
    pandas.Series, pandas.DataFrame or None
        A series {rid: flux} giving the import flux for each required import
        reaction and (optionally) the associated export fluxes. All exchange
        fluxes are oriented into the import reaction e.g. positive fluxes
        denote imports and negative fluxes exports. If `minimize_components`
        is a number larger 1 may return a DataFrame where each column is a
        minimal medium. Returns None if the minimization is infeasible
        (for instance if min_growth > maximum growth rate).

    Notes
    -----
    Due to numerical issues the `minimize_components` option will usually only
    minimized the number or "large" import fluxes. Specifically, the detection
    limit is given by `integrality_tolerance * max_bound` where `max_bound` is
    the largest bound on an import reaction. Thus, if you are interested
    in small import fluxes as well you may have to adjust the integrality
    tolerance at first with
    `model.solver.configuration.tolerances.integrality = 1e-7` for instance.
    However, this will be *very* slow for large models especially with GLPK.

    """
    LOGGER.info("calculating minimal medium for %s" % model.id)
    exchange_rxns = exchanges(model)
    if isinstance(open_exchanges, bool):
        open_bound = 1000
    else:
        open_bound = open_exchanges

    with model as mod:
        if open_exchanges:
            LOGGER.info("opening exchanges for %d imports" %
                        len(exchange_rxns))
            for rxn in exchange_rxns:
                rxn.bounds = (-open_bound, open_bound)
        LOGGER.info("applying growth rate constraints")
        obj_const = mod.problem.Constraint(
            mod.objective.expression, lb=min_growth,
            name="medium_obj_constraint")
        mod.add_cons_vars([obj_const])
        mod.solver.update()
        mod.objective = Zero
        LOGGER.info("adding new media objective")
        tol = mod.solver.configuration.tolerances.feasibility

        if minimize_components:
            add_mip_obj(mod)
            if isinstance(minimize_components, bool):
                minimize_components = 1
            seen = set()
            best = num_components = mod.slim_optimize()
            if mod.solver.status != OPTIMAL:
                LOGGER.warning("minimization of medium was infeasible")
                return None
            exclusion = mod.problem.Constraint(Zero, ub=0)
            mod.add_cons_vars([exclusion])
            mod.solver.update()
            media = []
            for i in range(minimize_components):
                LOGGER.info("finding alternative medium #%d" % (i + 1))
                vars = [mod.variables["ind_" + s] for s in seen]
                if len(seen) > 0:
                    exclusion.set_linear_coefficients(
                        dict.fromkeys(vars, 1))
                    exclusion.ub = best - 1
                num_components = mod.slim_optimize()
                if mod.solver.status != OPTIMAL or num_components > best:
                    break
                medium = _as_medium(exchange_rxns, tol, exports=exports)
                media.append(medium)
                seen.update(medium.index[medium > 0])
            if len(media) > 1:
                medium = pd.concat(media, axis=1).fillna(0.0)
            else:
                medium = media[0]
        else:
            add_linear_obj(mod)
            mod.slim_optimize()
            if mod.solver.status != OPTIMAL:
                LOGGER.warning("minimization of medium was infeasible")
                return None
            medium = _as_medium(exchange_rxns, tol, exports=exports)

    return medium
