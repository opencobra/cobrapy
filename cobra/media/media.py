"""Manages functions for growth media analysis and manipulation."""

from collections import Counter
from sympy.core.singleton import S
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
    for rxn in model.exchanges:
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
    boundary_rxns = model.exchanges
    M = max(np.max(np.abs(r.bounds)) for r in boundary_rxns)
    prob = model.problem
    coefs = {}
    to_add = []
    for rxn in boundary_rxns:
        export = len(rxn.reactants) == 1
        indicator = prob.Variable("ind_" + rxn.id, lb=0, ub=1, type="binary")
        if export:
            vrv = rxn.reverse_variable
            indicator_const = prob.Constraint(
                vrv - indicator * M, ub=0, name="ind_constraint_" + rxn.name)
        else:
            vfw = rxn.forward_variable
            indicator_const = prob.Constraint(
                vfw - indicator * M, ub=0, name="ind_constraint_" + rxn.name)
        to_add.extend([indicator, indicator_const])
        coefs[indicator] = 1
    model.add_cons_vars(to_add)
    model.solver.update()
    model.objective.set_linear_coefficients(coefs)
    model.objective.direction = "min"


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
    minimize_components : boolean
        Whether to minimize the number of components instead of the total
        import flux. Might be more intuitive if set to True but may also be
        slow to calculate for large communities.
    open_exchanges : boolean or number
        Whether to ignore currently set bounds and make all exchange reactions
        in the model possible. If set to a number all exchange reactions will
        be opened with (-number, number) as bounds.

    Returns
    -------
    pandas.Series or None
        A series {rid: flux} giving the import flux for each required import
        reaction and (optionally) the associated export fluxes. All exchange
        fluxes are oriented into the import reaction e.g. positive fluxes
        denote imports and negative fluxes exports. Returns None if the
        minimization is infeasible.

    """
    LOGGER.info("calculating minimal medium for %s" % model.id)
    boundary_rxns = exchanges(model)
    if isinstance(open_exchanges, bool):
        open_bound = 1000
    else:
        open_bound = open_exchanges

    with model as mod:
        if open_exchanges:
            LOGGER.info("opening exchanges for %d imports" %
                        len(boundary_rxns))
            for rxn in boundary_rxns:
                rxn.bounds = (-open_bound, open_bound)
        LOGGER.info("applying growth rate constraints")
        obj_const = mod.problem.Constraint(
            mod.objective.expression, lb=min_growth,
            name="medium_obj_constraint")
        mod.add_cons_vars([obj_const])
        mod.solver.update()
        mod.objective = S.Zero
        LOGGER.info("adding new media objective")
        if minimize_components:
            add_mip_obj(mod)
        else:
            add_linear_obj(mod)
        mod.solver.optimize()
        if mod.solver.status != OPTIMAL:
            LOGGER.warning("minimization of medium was infeasible")
            return None

        LOGGER.info("formatting medium")
        medium = pd.Series()
        tol = model.solver.configuration.tolerances.feasibility
        for rxn in boundary_rxns:
            export = len(rxn.reactants) == 1
            flux = rxn.flux
            if abs(flux) < tol:
                continue
            if export:
                medium[rxn.id] = -flux
            elif not export:
                medium[rxn.id] = flux
        if not exports:
            medium = medium[medium > 0]

    return medium
