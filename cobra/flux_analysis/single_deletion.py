from warnings import warn

from six import string_types, iteritems

from ..manipulation import delete_model_genes, undelete_model_genes
from ..manipulation.delete import find_gene_knockout_reactions
from ..solvers import solver_dict, get_solver_name
try:
    import scipy
except ImportError:
    moma = None
else:
    from . import moma


def single_deletion(cobra_model, element_list=None,
                    element_type='gene', **kwargs):
    """Wrapper for single_gene_deletion and single_reaction_deletion

    .. deprecated :: 0.4
        Use single_reaction_deletion and single_gene_deletion
    """
    warn("deprecated - use single_reaction_deletion and single_gene_deletion")
    if element_type == "reaction":
        return single_reaction_deletion(cobra_model, element_list, **kwargs)
    elif element_type == "gene":
        return single_gene_deletion(cobra_model, element_list, **kwargs)
    else:
        raise Exception("unknown element type")


def single_reaction_deletion(cobra_model, reaction_list=None, solver=None,
                             method="fba", **solver_args):
    """sequentially knocks out each reaction in a model

    reaction_list: list of reaction_ids or cobra.Reaction

    method: "fba" or "moma"

    returns ({reaction_id: growth_rate}, {reaction_id: status})"""
    if reaction_list is None:
        reaction_list = cobra_model.reactions
    else:
        reaction_list = [cobra_model.reactions.get_by_id(i)
                         if isinstance(i, string_types) else i
                         for i in reaction_list]
    if method == "fba":
        return single_reaction_deletion_fba(cobra_model, reaction_list,
                                            solver=solver, **solver_args)
    elif method == "moma":
        return single_reaction_deletion_moma(cobra_model, reaction_list,
                                             solver=solver, **solver_args)
    else:
        raise ValueError("Unknown deletion method '%s'" % method)


def single_reaction_deletion_fba(cobra_model, reaction_list, solver=None,
                                 **solver_args):
    """sequentially knocks out each reaction in a model using FBA

    reaction_list: list of reaction_ids or cobra.Reaction

    method: "fba" or "moma"

    returns ({reaction_id: growth_rate}, {reaction_id: status})"""

    solver = solver_dict[get_solver_name() if solver is None else solver]
    lp = solver.create_problem(cobra_model)

    growth_rate_dict = {}
    status_dict = {}
    for reaction in reaction_list:
        old_bounds = (reaction.lower_bound, reaction.upper_bound)
        index = cobra_model.reactions.index(reaction)
        solver.change_variable_bounds(lp, index, 0., 0.)
        solver.solve_problem(lp, **solver_args)
        # get the status and growth rate
        status = solver.get_status(lp)
        status_dict[reaction.id] = status
        growth_rate_dict[reaction.id] = solver.get_objective_value(lp) \
            if status == "optimal" else 0.
        # reset the problem
        solver.change_variable_bounds(lp, index, old_bounds[0], old_bounds[1])
    return (growth_rate_dict, status_dict)


def single_reaction_deletion_moma(cobra_model, reaction_list, solver=None,
                                  **solver_args):
    """sequentially knocks out each reaction in a model using MOMA

    reaction_list: list of reaction_ids or cobra.Reaction


    returns ({reaction_id: growth_rate}, {reaction_id: status})"""
    # The same function can not be used because MOMA can not re-use the
    # same LP object. Problem re-use leads to incorrect solutions.
    if moma is None:
        raise RuntimeError("scipy required for moma")
    solver = solver_dict[solver if solver else get_solver_name(qp=True)]
    moma_model, moma_objective = moma.create_euclidian_moma_model(cobra_model)

    growth_rate_dict = {}
    status_dict = {}
    for reaction in reaction_list:
        index = cobra_model.reactions.index(reaction)
        solution = moma.moma_knockout(moma_model, moma_objective, (index,),
                                      solver=solver, **solver_args)
        status_dict[reaction.id] = solution.status
        growth_rate_dict[reaction.id] = solution.f
    return (growth_rate_dict, status_dict)


def single_gene_deletion(cobra_model, gene_list=None, solver=None,
                         method="fba", **solver_args):
    """sequentially knocks out each gene in a model

    gene_list: list of gene_ids or cobra.Gene

    method: "fba" or "moma"

    returns ({gene_id: growth_rate}, {gene_id: status})"""
    if gene_list is None:
        gene_list = cobra_model.genes
    else:
        gene_list = [cobra_model.genes.get_by_id(i)
                     if isinstance(i, string_types) else i for i in gene_list]

    if method == "fba":
        return single_gene_deletion_fba(cobra_model, gene_list,
                                        solver=solver, **solver_args)
    elif method == "moma":
        return single_gene_deletion_moma(cobra_model, gene_list,
                                         solver=solver, **solver_args)
    else:
        raise ValueError("Unknown deletion method '%s'" % method)


def single_gene_deletion_fba(cobra_model, gene_list, solver=None,
                             **solver_args):

    solver = solver_dict[get_solver_name() if solver is None else solver]
    lp = solver.create_problem(cobra_model)

    growth_rate_dict = {}
    status_dict = {}
    for gene in gene_list:
        old_bounds = {}
        for reaction in find_gene_knockout_reactions(cobra_model, [gene]):
            index = cobra_model.reactions.index(reaction)
            old_bounds[index] = (reaction.lower_bound, reaction.upper_bound)
            solver.change_variable_bounds(lp, index, 0., 0.)
        solver.solve_problem(lp, **solver_args)
        # get the status and growth rate
        status = solver.get_status(lp)
        status_dict[gene.id] = status
        growth_rate = solver.get_objective_value(lp) \
            if status == "optimal" else 0.
        growth_rate_dict[gene.id] = growth_rate
        # reset the problem
        for index, bounds in iteritems(old_bounds):
            solver.change_variable_bounds(lp, index, bounds[0], bounds[1])
    return (growth_rate_dict, status_dict)


def single_gene_deletion_moma(cobra_model, gene_list, solver=None,
                              **solver_args):
    if moma is None:
        raise RuntimeError("scipy required for moma")
    solver = solver if solver else get_solver_name(qp=True)
    moma_model, moma_objective = moma.create_euclidian_moma_model(cobra_model)

    growth_rate_dict = {}
    status_dict = {}
    for gene in gene_list:
        delete_model_genes(moma_model, [gene.id])
        solution = moma.solve_moma_model(moma_model, moma_objective,
                                         solver=solver, **solver_args)
        status_dict[gene.id] = solution.status
        growth_rate_dict[gene.id] = solution.f
        undelete_model_genes(moma_model)
    return (growth_rate_dict, status_dict)
