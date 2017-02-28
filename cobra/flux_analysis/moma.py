# -*- coding: utf-8 -*-

from __future__ import absolute_import

from scipy.sparse import dok_matrix
from sympy.core.singleton import S

import cobra.util.solver as sutil
from cobra.solvers import get_solver_name, solver_dict


def add_moma(model):
    """Add constraints and objective representing a MOMA
    (minimization of metabolic adjustment) model.

    Parameters:
    -----------

    model : a cobra model

    Returns:
    --------
    Nothing.

    Notes
    -----

    In the original MOMA specification one looks for the flux distribution
    of the deletion (v^d) closest to the fluxes without the deletion (v).
    In math this means:

    minimize \sum_i (v^d_i - v_i)^2
    s.t. Sv^d = 0
         lb_i <= v^d_i <= ub_i

    Here, we use a variable transformation v^t := v^d_i - v_i. Substituting
    and using the fact that Sv = 0 gives:

    minimize \sum (v^t_i)^2
    s.t. Sv^d = 0
         v^t = v^d_i - v_i
         lb_i <= v^d_i <= ub_i

    So basically we just re-center the flux space at the old solution and than
    find the flux distribution closest to the new zero (center). This is the
    same strategy as used in cameo.

    The former objective function is saved in the optlang solver interface as
    "moma_old_objective" and this can be used to immediately extract the value
    of the former objective after MOMA optimization.
    """
    if 'moma_old_objective' in model.solver.variables:
        raise ValueError('model is already adjusted for MOMA')

    # Fall back to default QP solver if current one has no QP capability
    model.solver = sutil.choose_solver(model, qp=True)[1]

    solution = model.optimize()
    prob = model.problem
    v = prob.Variable("moma_old_objective")
    c = prob.Constraint(model.solver.objective.expression - v,
                        lb=0.0, ub=0.0, name="moma_old_objective_constraint")
    to_add = [v, c]
    new_obj = S.Zero
    for r in model.reactions:
        flux = solution.fluxes[r.id]
        dist = prob.Variable("moma_dist_" + r.id)
        const = prob.Constraint(r.flux_expression - dist, lb=flux, ub=flux,
                                name="moma_constraint_" + r.id)
        to_add.extend([dist, const])
        new_obj += dist**2
    model.add_cons_vars(to_add)
    model.objective = prob.Objective(new_obj, direction='min')


def create_euclidian_moma_model(cobra_model, wt_model=None, **solver_args):
    # make the wild type copy if none was supplied
    if wt_model is None:
        wt_model = cobra_model.copy()
    else:
        wt_model = wt_model.copy()
        # ensure single objective
        wt_obj = sutil.linear_reaction_coefficients(wt_model)
        if len(wt_obj) != 1:
            raise ValueError("wt_model must have exactly 1 objective, %d found"
                             % len(wt_obj))

    obj = sutil.linear_reaction_coefficients(wt_model)
    if len(obj) == 1:
        objective_id = list(obj)[0].id
    else:
        raise ValueError("model must have exactly 1 objective, %d found" %
                         len(obj))

    wt_model.optimize(**solver_args)
    for reaction in wt_model.reactions:
        # we don't want delete_model_gene to remove the wt reaction as well
        reaction.gene_reaction_rule = ''
        if reaction.objective_coefficient != 0:
            reaction.objective_coefficient = 0
            reaction.upper_bound = reaction.lower_bound = reaction.x
        reaction.id = "MOMA_wt_" + reaction.id
    for metabolite in wt_model.metabolites:
        metabolite.id = "MOMA_wt_" + metabolite.id
    wt_model.repair()

    # make the moma model by combining both
    moma_model = cobra_model.copy()
    for reaction in moma_model.reactions:
        reaction.objective_coefficient = 0
    moma_model.add_reactions(wt_model.reactions)
    return moma_model, objective_id


def create_euclidian_distance_objective(n_moma_reactions):
    """returns a matrix which will minimze the euclidian distance

    This matrix has the structure
    [ I  -I]
    [-I   I]
    where I is the identity matrix the same size as the number of
    reactions in the original model.

    n_moma_reactions: int
        This is the number of reactions in the MOMA model, which should
        be twice the number of reactions in the original model"""
    if n_moma_reactions % 2 != 0:
        raise ValueError("must be even")
    n_reactions = n_moma_reactions // 2
    Q = dok_matrix((n_reactions * 2, n_reactions * 2))
    for i in range(2 * n_reactions):
        Q[i, i] = 1
    for i in range(n_reactions):
        Q[i, n_reactions + i] = -1
        Q[n_reactions + i, i] = -1
    return Q


def create_euclidian_distance_lp(moma_model, solver):
    Q = create_euclidian_distance_objective(len(moma_model.reactions))
    lp = solver.create_problem(moma_model, objective_sense="minimize",
                               quadratic_component=Q)
    return lp


def solve_moma_model(moma_model, objective_id, solver=None, **solver_args):
    solver = solver_dict[solver if solver and isinstance(solver, str)
                         else get_solver_name(qp=True)]
    lp = create_euclidian_distance_lp(moma_model, solver=solver)
    solver.solve_problem(lp, **solver_args)
    solution = solver.format_solution(lp, moma_model)
    solution.f = 0. if solution.x_dict is None \
        else solution.x_dict[objective_id]
    return solution


def moma(wt_model, mutant_model, solver=None, **solver_args):
    if "norm_type" in solver_args:
        print("only euclidian norm type supported for moma")
        solver_args.pop("norm_type")
    moma_model, objective_id = create_euclidian_moma_model(mutant_model,
                                                           wt_model)
    return solve_moma_model(moma_model, objective_id,
                            solver=solver, **solver_args)


def moma_knockout(moma_model, moma_objective, reaction_indexes, **moma_args):
    """computes result of reaction_knockouts using moma"""
    n = len(moma_model.reactions) // 2
    # knock out the reaction
    for i in reaction_indexes:
        mutant_reaction = moma_model.reactions[i]
        mutant_reaction.lower_bound, mutant_reaction.upper_bound = (0., 0.)
    result = solve_moma_model(moma_model, moma_objective, **moma_args)
    # reset the knockouts
    for i in reaction_indexes:
        mutant_reaction = moma_model.reactions[i]
        wt_reaction = moma_model.reactions[n + i]
        mutant_reaction.lower_bound = wt_reaction.lower_bound
        mutant_reaction.upper_bound = wt_reaction.upper_bound
    return result
