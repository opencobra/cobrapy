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

from collections import Counter
import logging

LOGGER = logging.getLogger(__name__)

excludes = {"demand": ["SN_", "SK_", "sink", "EX_", "exchange"],
            "exchange": ["demand", "DM_", "biosynthesis", "transcription",
                         "replication", "SN_", "SK_", "sink"],
            "sink": ["demand", "DM_", "biosynthesis", "transcription",
                     "replication", "EX_", "exchange"]}
"""A list of sub-strings in reaction IDs that usually indicate
that the reaction is *not* a reaction of the specified type."""

sbo_terms = {"demand": "SBO:0000628",
             "exchange": "SBO:0000627",
             "sink": "SBO:0000632",
             "biomass": "SBO:0000629",
             "pseudoreaction": "SBO:0000631"}
"""SBO term identifiers for various boundary types."""


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
    if not model.boundary:
        LOGGER.error("The heuristic for discovering an external compartment "
                     "relies on boundary reactions. Yet, there are no "
                     "boundary reactions in this model.")
        raise RuntimeError(
            "The external compartment cannot be identified. "
            "The heuristic for discovering an external compartment "
            "relies on boundary reactions. Yet, there are no "
            "boundary reactions in this model.")
    counts = Counter(tuple(r.compartments)[0] for r in model.boundary)
    most = counts.most_common(1)[0][0]
    if "e" in model.compartments:
        if most == "e":
            return "e"
        else:
            LOGGER.warning("There is an `e` compartment but it does not look "
                           "like it is the actual external compartment.")
        return most
    return most


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
    sbo_term = reaction.annotation.get("SBO", "").upper()
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
