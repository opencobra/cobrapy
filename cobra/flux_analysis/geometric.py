# -*- coding: utf8 -*-

"Contains functions to perform geometric FBA."

from __future__ import absolute_import, division

from optlang.symbolics import Zero

from cobra.flux_analysis.parsimonious import add_pfba
from cobra.flux_analysis.variability import flux_variability_analysis
from cobra.core import get_solution


def geometric_fba(model, epsilon=1E-06, max_tries=200):
    """Perform geometric FBA to obtain a unique, centered flux distribution.

    Geometric FBA [1]_ formulates the problem as a polyhedron and
    then solves it by bounding the convex hull of the polyhedron.
    The bounding forms a box around the convex hull which reduces
    with every iteration and extracts a unique solution in this way.

    Parameters
    ----------
    model: cobra.Model
        The model to perform geometric FBA on.
    epsilon: float, optional
        The convergence tolerance of the model (default 1E-06).
    max_tries: int, optional
        Maximum number of iterations (default 200).

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
        # iteration parameters
        delta = 1.0  # initialize at 1.0 to enter while loop
        count = 2  # iteration #1 happens out of the loop

        # vars and consts storage variables
        consts = []
        obj_vars = []
        updating_vars_cons = []

        # first iteration
        prob = model.problem
        add_pfba(model)  # minimizes the solution space to convex hull
        model.optimize()
        fva_sol = flux_variability_analysis(model)
        mean_flux = (fva_sol["maximum"] + fva_sol["minimum"]).abs() / 2

        # set gFBA constraints
        for rxn in model.reactions:
            var = prob.Variable("geometric_fba_" + rxn.id,
                                lb=0,
                                ub=mean_flux[rxn.id])
            upper_const = prob.Constraint(rxn.flux_expression - var,
                                          ub=mean_flux[rxn.id],
                                          name="geometric_fba_upper_const_" +
                                          rxn.id)
            lower_const = prob.Constraint(rxn.flux_expression + var,
                                          lb=fva_sol.at[rxn.id, "minimum"],
                                          name="geometric_fba_lower_const_" +
                                          rxn.id)
            updating_vars_cons.append((rxn.id, var, upper_const, lower_const))
            consts.extend([var, upper_const, lower_const])
            obj_vars.append(var)
        model.add_cons_vars(consts)

        # minimize distance between flux and centre
        model.objective = prob.Objective(Zero, sloppy=True, direction="min")
        model.objective.set_linear_coefficients({v: 1.0 for v in obj_vars})
        model.optimize()

        # further iterations
        while delta > epsilon and count <= max_tries:
            fva_sol = flux_variability_analysis(model)
            mean_flux = (fva_sol["maximum"] + fva_sol["minimum"]).abs() / 2

            for rxn_id, var, u_c, l_c in updating_vars_cons:
                var.ub = mean_flux[rxn_id]
                u_c.ub = mean_flux[rxn_id]
                l_c.lb = fva_sol.at[rxn_id, "minimum"]
            model.optimize()
            delta = (fva_sol["maximum"] - fva_sol["minimum"]).max()

            count += 1
            if count == max_tries:
                raise RuntimeError(
                    "The iterations have exceeded the maximum value of {}. "
                    "This is probably due to the increased complexity of the "
                    "model and can lead to inaccurate results. Please set a "
                    "different convergence tolerance and/or increase the "
                    "maximum iterations".format(max_tries)
                )

    return get_solution(model)
