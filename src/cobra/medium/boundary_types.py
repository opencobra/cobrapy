"""Provide functions to identify the type of boundary reactions.

This module uses various heuristics to decide whether a boundary reaction
is an exchange, demand or sink reaction. It mostly orientates on the
following paper:

Thiele, I., & Palsson, B. Ø. (2010, January). A protocol for
generating a high-quality genome-scale metabolic reconstruction.
Nature protocols. Nature Publishing Group.
http://doi.org/10.1038/nprot.2009.203

"""

import logging
from typing import TYPE_CHECKING, List, Optional

import pandas as pd

from .annotations import compartment_shortlist, excludes, sbo_terms


if TYPE_CHECKING:
    from cobra import Model, Reaction


logger = logging.getLogger(__name__)


def find_external_compartment(model: "Model") -> str:
    """Find the external compartment in the model.

    Uses a simple heuristic where the external compartment should be the
    one with the most exchange reactions.

    Parameters
    ----------
    model : cobra.Model
        The cobra model whose external compartments are to be identified.

    Returns
    -------
    str
        The putative external compartment.

    Raises
    ------
    RuntimeError
        If several compartments are similar and thus difficult to identify,
        or, recognized names usually used for external compartment are
        absent.

    """
    if model.boundary:
        counts = pd.Series(tuple(r.compartments)[0] for r in model.boundary)
        most = counts.value_counts()
        most = most.index[most == most.max()].to_series()
    else:
        most = None
    like_external = compartment_shortlist["e"] + ["e"]
    matches = pd.Series(
        [co in like_external for co in model.compartments],
        dtype=bool,
        index=model.compartments,
    )

    if matches.sum() == 1:
        compartment = matches.index[matches][0]
        logger.info(
            f"Compartment `{compartment}` sounds like an external compartment. "
            "Using this one without counting boundary reactions."
        )
        return compartment
    elif most is not None and matches.sum() > 1 and matches[most].sum() == 1:
        compartment = most[matches[most]][0]
        logger.warning(
            "There are several compartments that look like an "
            f"external compartment but `{compartment}` has the most boundary "
            "reactions, so using that as the external compartment."
        )
        return compartment
    elif matches.sum() > 1:
        raise RuntimeError(
            "There are several compartments that look "
            "like external compartments but we can't tell "
            "which one to use. Consider renaming your "
            "compartments please."
        )

    if most is not None:
        logger.warning(
            "Could not identify an external compartment by name and "
            "choosing one with the most boundary reactions. That "
            "might be complete nonsense or change suddenly. "
            "Consider renaming your compartments using "
            "`Model.compartments` to fix this."
        )
        return most[0]

    # No info in the model, so give up
    raise RuntimeError(
        "The heuristic for discovering an external compartment "
        "relies on names and boundary reactions. Yet, there "
        "are neither compartments with recognized names nor "
        "boundary reactions in the model."
    )


def is_boundary_type(
    reaction: "Reaction", boundary_type: str, external_compartment: str
) -> bool:
    """Check whether a reaction is an exchange reaction.

    Parameters
    ----------
    reaction : cobra.Reaction
        The reaction to check.
    boundary_type : {"exchange", "demand", "sink"}
        Boundary type to check for.
    external_compartment : str
        The ID for the external compartment.

    Returns
    -------
    bool
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

    return (
        reaction.boundary
        and not any(ex in reaction.id for ex in excludes[boundary_type])
        and correct_compartment
        and rev_type
    )


def find_boundary_types(
    model: "Model", boundary_type: str, external_compartment: Optional[str] = None
) -> List["Reaction"]:
    """Find specific boundary reactions.

    Parameters
    ----------
    model : cobra.Model
        The cobra model whose boundary reactions are to be found.
    boundary_type : {"exchange", "demand", "sink"}
        Boundary type to check for.
    external_compartment : str, optional
        The ID for the external compartment. If None, it will be detected
        automatically (default None).

    Returns
    -------
    list of cobra.Reaction or an empty list
        A list of likely boundary reactions of a user defined type.

    """
    if not model.boundary:
        logger.warning(
            "There are no boundary reactions in this model. "
            "Therefore, specific types of boundary reactions such "
            "as 'exchanges', 'demands' or 'sinks' cannot be "
            "identified."
        )
        return []
    if external_compartment is None:
        external_compartment = find_external_compartment(model)
    return model.reactions.query(
        lambda r: is_boundary_type(r, boundary_type, external_compartment)
    )
