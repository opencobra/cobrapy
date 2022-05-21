"""Provide function for annotating demand and exchange reactions."""

from typing import TYPE_CHECKING


if TYPE_CHECKING:
    from cobra import Model


def add_SBO(model: "Model") -> None:
    """Add SBO terms for demands and exchanges.

    This works for models which follow the standard convention for
    constructing and naming these reactions.

    The reaction should only contain the single metabolite being exchanged,
    and the id should be EX_<met_id> or DM_<met_id> .

    Parameters
    ----------
    model: cobra.Model
        The model whose demand and exchange reactions need to be annotated.

    """
    #??? Should this be done with boundary?
    for r in model.reactions:
        # don't annotate already annotated reactions
        if len(r.annotation.get("sbo")) != 0 and r.annotation.sbo:
            continue
        # only doing exchanges
        if len(r.metabolites) != 1:
            continue
        met_id = list(r._metabolites)[0].id
        if r.id.startswith("EX_") and r.id == "EX_" + met_id:
            r.annotation["sbo"] = ["SBO:0000627"]
        elif r.id.startswith("DM_") and r.id == "DM_" + met_id:
            r.annotation["sbo"] = ["SBO:0000628"]
        elif r.id.startswith("sink_") and r.id == "Sink_" + met_id:
            r.annotation["sbo"] = ["SBO:0000632"]
