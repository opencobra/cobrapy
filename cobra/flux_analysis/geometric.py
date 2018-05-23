# -*- coding: utf8 -*-

"Contains functions to perform geometric FBA."

from __future__ import absolute_import

from optlang.symbolics import Zero

from cobra.flux_analysis.parsimonious import add_pfba
from cobra.flux_analysis.variability import flux_variability_analysis
from cobra.core import get_solution


def geometric_fba(model, epsilon=1E-06):
    """Perform geometric FBA to obtain a unique
    centered flux.

    Geometric FBA [1] formulates the problem as a polyhedron and
    then solves it by bounding the convex hull of the polyhedron.
    The bounding forms a box around the convex hull which reduces
    with every iteration and extracts a unique solution in this way.

    Parameters
    ----------
    model: cobra.Model
        The model to perform geometric FBA on.
    epsilon: float
        The convergence tolerance of the model.
        (default 10e-6).

    Returns
    -------
    cobra.Solution
        The solution object containing all the constraints required
        for geometric FBA.

    References
    ----------
    .. [1] Smallbone, Kieran & Simeonidis, Vangelis. (2009).
    Flux balance analysis: A geometric perspective.
    Journal of theoretical biology.258. 311-5. 10.1016/j.jtbi.2009.01.027.
    """

    with model:
        prob = model.problem
        add_pfba(model)  # minimizes the solution space to convex hull
        delta = 1.0  # initialize at 1.0 to enter while loop
        # first iteration
        sol = model.optimize()
        fva_sol = flux_variability_analysis(model)
        v_n = 0.5 * abs(fva_sol['maximum'] + fva_sol['minimum'])
        consts = []
        obj_vars = []
        # set gFBA constraints
        for rxn in model.reactions:
            var = prob.Variable("geometric_fba_" + rxn.id,
                                lb=0,
                                ub=v_n[rxn.id])
            upper_const = prob.Constraint(rxn.flux_expression - var,
                                          ub=v_n[rxn.id],
                                          name="geometric_fba_upper_const_" +
                                          rxn.id)
            lower_const = prob.Constraint(rxn.flux_expression + var,
                                          lb=fva_sol['minimum'][rxn.id],
                                          name="geometric_fba_lower_const_" +
                                          rxn.id)
            consts.extend([var, upper_const, lower_const])
            obj_vars.append(var)
        model.add_cons_vars(consts)
        # minimize distance between flux and centre
        model.objective = prob.Objective(Zero, sloppy=True, direction="min")
        model.objective.set_linear_coefficients({v: 1.0 for v in obj_vars})
        model.optimize()
        while delta > epsilon:
            fva_sol = (flux_variability_analysis(model))
            v_n = 0.5 * abs(fva_sol['maximum'] + fva_sol['minimum'])
            for rxn in model.reactions:
                model.variables["geometric_fba_" + rxn.id].ub = v_n[rxn.id]
                model.constraints["geometric_fba_upper_const_" + rxn.id].ub = \
                    v_n[rxn.id]
                model.constraints["geometric_fba_lower_const_" + rxn.id].lb = \
                    fva_sol['minimum'][rxn.id]
            model.optimize()
            delta = max(fva_sol['maximum'] - fva_sol['minimum'])

    return get_solution(model)
