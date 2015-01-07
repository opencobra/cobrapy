from ..manipulation import modify


def optimize_minimal_flux(model, already_irreversible=False,
                          **optimize_kwargs):
    """Perform basic pFBA (parsimonius FBA) and minimize total flux.

    The function attempts to act as a drop-in replacement for optimize. It
    will make the reaction reversible and perform an optimization, then
    force the objective value to remain the same and minimize the total
    flux. Finally, it will convert the reaction back to the irreversible
    form it was in before. See http://dx.doi.org/10.1038/msb.2010.47

    model : :class:`~cobra.core.Model` object

    already_irreversible : bool, optional
        By default, the model is converted to an irreversible one.
        However, if the model is already irreversible, this step can be
        skipped.

    """
    if "new_objective" in optimize_kwargs:
        raise ValueError("Use objective coefficients, not new_objective")
    if not already_irreversible:
        modify.convert_to_irreversible(model)
    model.optimize(**optimize_kwargs)
    # if the problem is infeasible
    if model.solution.f is None:
        raise Exception("model could not be solved")
    old_f = model.solution.f
    old_objective_coefficients = {}
    old_lower_bounds = {}
    old_upper_bounds = {}
    for reaction in model.reactions:
        # if the reaction has a nonzero objective coefficient, then
        # the same flux should be maintained through that reaction
        if reaction.objective_coefficient != 0:
            old_objective_coefficients[reaction] = \
                reaction.objective_coefficient
            old_lower_bounds[reaction] = reaction.lower_bound
            old_upper_bounds[reaction] = reaction.upper_bound
            x = model.solution.x_dict[reaction.id]
            reaction.lower_bound = x
            reaction.upper_bound = x
            reaction.objective_coefficient = 0
        else:
            reaction.objective_coefficient = 1
    # set to minimize flux
    optimize_kwargs["objective_sense"] = "minimize"
    model.optimize(**optimize_kwargs)
    # make the model back the way it was
    for reaction in model.reactions:
        if reaction in old_objective_coefficients:
            reaction.objective_coefficient = \
                old_objective_coefficients[reaction]
            reaction.lower_bound = old_lower_bounds[reaction]
            reaction.upper_bound = old_upper_bounds[reaction]
        else:
            reaction.objective_coefficient = 0
    # if the minimization problem was successful
    if model.solution.f is not None:
        model.solution.f = old_f
    if not already_irreversible:
        modify.revert_to_reversible(model)
    return model.solution
