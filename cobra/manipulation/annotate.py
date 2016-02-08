from six import iteritems


def add_SBO(model):
    """adds SBO terms for demands and exchanges

    This works for models which follow the standard convention for
    constructing and naming these reactions.

    The reaction should only contain the single metabolite being exchanged,
    and the id should be EX_metid or DM_metid
    """
    for r in model.reactions:
        # don't annotate already annotated reactions
        if r.annotation.get("SBO"):
            continue
        # only doing exchanges
        if len(r.metabolites) != 1:
            continue
        met_id = list(r._metabolites)[0].id
        if r.id.startswith("EX_") and r.id == "EX_" + met_id:
            r.annotation["SBO"] = "SBO:0000627"
        elif r.id.startswith("DM_") and r.id == "DM_" + met_id:
            r.annotation["SBO"] = "SBO:0000628"
