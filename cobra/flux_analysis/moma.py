from scipy.sparse import dok_matrix

from ..solvers import get_solver_name, solver_dict


def create_euclidian_moma_model(cobra_model, wt_model=None, **solver_args):
    # make the wild type copy if none was supplied
    if wt_model is None:
        wt_model = cobra_model.copy()
    else:
        wt_model = wt_model.copy()
        # ensure single objective
        wt_obj = wt_model.reactions.query(lambda x: x > 0,
                                          "objective_coefficient")
        if len(wt_obj) != 1:
            raise ValueError("wt_model must have exactly 1 objective, %d found"
                             % len(wt_obj))

    obj = cobra_model.reactions.query(lambda x: x > 0, "objective_coefficient")
    if len(obj) == 1:
        objective_id = obj[0].id
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
    moma_model.solution = solution
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
