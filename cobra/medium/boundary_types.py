# -*- coding: utf-8 -*-

"""Contains function to identify the type of boundary reactions.

This module uses various heuristics to decide whether a boundary reaction
is an exchange, demand or sink reaction. It mostly orientates on the following
paper:

Thiele, I., & Palsson, B. Ø. (2010, January). A protocol for
generating a high-quality genome-scale metabolic reconstruction.
Nature protocols. Nature Publishing Group.
http://doi.org/10.1038/nprot.2009.203
"""

import logging
from collections import Counter

import pandas as pd

from cobra.medium.annotations import compartment_shortlist, excludes, sbo_terms


LOGGER = logging.getLogger(__name__)


def find_external_compartment(model):
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
    if model.boundary:
        counts = pd.Series(tuple(r.compartments)[0] for r in model.boundary)
        most = counts.value_counts()
        most = most.index[most == most.max()].to_series()
    else:
        most = None
    like_external = compartment_shortlist["e"] + ["e"]
    matches = pd.Series([co in like_external for co in model.compartments],
                        index=model.compartments)

    if matches.sum() == 1:
        compartment = matches.index[matches][0]
        LOGGER.info("Compartment `%s` sounds like an external compartment. "
                    "Using this one without counting boundary reactions" %
                    compartment)
        return compartment
    elif most is not None and matches.sum() > 1 and matches[most].sum() == 1:
        compartment = most[matches[most]][0]
        LOGGER.warning("There are several compartments that look like an "
                       "external compartment but `%s` has the most boundary "
                       "reactions, so using that as the external "
                       "compartment." % compartment)
        return compartment
    elif matches.sum() > 1:
        raise RuntimeError("There are several compartments (%s) that look "
                           "like external compartments but we can't tell "
                           "which one to use. Consider renaming your "
                           "compartments please.")

    if most is not None:
        return most[0]
        LOGGER.warning("Could not identify an external compartment by name and"
                       " choosing one with the most boundary reactions. That "
                       "might be complete nonsense or change suddenly. "
                       "Consider renaming your compartments using "
                       "`Model.compartments` to fix this.")
    # No info in the model, so give up
    raise RuntimeError("The heuristic for discovering an external compartment "
                       "relies on names and boundary reactions. Yet, there "
                       "are neither compartments with recognized names nor "
                       "boundary reactions in the model.")


def is_boundary_type(reaction, boundary_type, external_compartment):
    """Check whether a reaction is an exchange reaction.

    Arguments
    ---------
    reaction : cobra.Reaction
        The reaction to check.
    boundary_type : str
        What boundary type to check for. Must be one of
        "exchange", "demand", or "sink".
    external_compartment : str
        The id for the external compartment.

    Returns
    -------
    boolean
        Whether the reaction looks like the requested type. Might be based
        on a heuristic.
    """
    # Check if the reaction has an annotation. Annotations dominate everything.
    sbo_term = reaction.annotation.get("sbo", "")
    if isinstance(sbo_term, list):
        sbo_term = sbo_term[0]
    sbo_term = sbo_term.upper()

    if sbo_term == sbo_terms[boundary_type]:
        return True
    if sbo_term in [sbo_terms[k] for k in sbo_terms if k != boundary_type]:
        return False

    # Check if the reaction is in the correct compartment (exterior or inside)
    correct_compartment = external_compartment in reaction.compartments
    if boundary_type != "exchange":
        correct_compartment = not correct_compartment

    # Check if the reaction has the correct reversibility
    rev_type = True
    if boundary_type == "demand":
        rev_type = not reaction.reversibility
    elif boundary_type == "sink":
        rev_type = reaction.reversibility

    return (reaction.boundary and not
            any(ex in reaction.id for ex in excludes[boundary_type]) and
            correct_compartment and rev_type)


def find_boundary_types(model, boundary_type, external_compartment=None):
    """Find specific boundary reactions.

    Arguments
    ---------
    model : cobra.Model
        A cobra model.
    boundary_type : str
        What boundary type to check for. Must be one of
        "exchange", "demand", or "sink".
    external_compartment : str or None
        The id for the external compartment. If None it will be detected
        automatically.

    Returns
    -------
    list of cobra.reaction
        A list of likely boundary reactions of a user defined type.
    """
    if not model.boundary:
        LOGGER.warning("There are no boundary reactions in this model. "
                       "Therefore specific types of boundary reactions such "
                       "as 'exchanges', 'demands' or 'sinks' cannot be "
                       "identified.")
        return []
    if external_compartment is None:
        external_compartment = find_external_compartment(model)
    return model.reactions.query(
        lambda r: is_boundary_type(r, boundary_type, external_compartment))
