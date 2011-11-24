import cobra


def optimize_minimal_flux(model, already_irreversible=False,
        **optimize_kwargs):
    """Perform basic pFBA (parsimonius FBA) and minimize flux
    through the network. The function attempts to act as a drop-in
    replacement for optimize."""
    if "new_objective" in optimize_kwargs:
        raise ValueError("Not implemented yet, use objective coefficients")
    if not already_irreversible:
        cobra.manipulation.modify.convert_to_irreversible(model)
    model.optimize(**optimize_kwargs)
    # if the problem is infeasible
    if model.solution.f is None:
        return
    old_f = model.solution.f
    old_objective_coefficients = {}
    old_lower_bounds = {}
    old_upper_bounds = {}
    for reaction in model.reactions:
        if reaction.objective_coefficient != 0:
            old_objective_coefficients[reaction] = \
                reaction.objective_coefficient
            old_lower_bounds[reaction] = reaction.lower_bound
            old_upper_bounds[reaction] = reaction.upper_bound
            reaction.lower_bound = reaction.x
            reaction.upper_bound = reaction.x
        reaction.objective_coefficient = -1
    model.optimize(**optimize_kwargs)
    # make the model back the way it was
    for reaction in old_objective_coefficients:
        reaction.objective_coefficient = old_objective_coefficients[reaction]
        reaction.lower_bound = old_lower_bounds[reaction]
        reaction.upper_bound = old_upper_bounds[reaction]
    model.solution.f = old_f
    cobra.manipulation.modify.convert_back_to_reversible(model)


if __name__ == "__main__":
    import cobra.test
    import cobra.oven.legacyIO
    model = cobra.oven.legacyIO.load_pickle(cobra.test.ecoli_pickle)
    optimize_minimal_flux(model)
