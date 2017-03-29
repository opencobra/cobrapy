# -*- coding: utf-8 -*-

from __future__ import absolute_import

import pandas
from six import iteritems, string_types
from optlang.interface import OPTIMAL

import cobra.solvers as legacy_solvers
import cobra.util.solver as solvers
from cobra.manipulation import delete_model_genes, undelete_model_genes
from cobra.manipulation.delete import find_gene_knockout_reactions

# this can be removed after deprecation of the old solver interface
# since the optlang vrsion requires neither numpy nor scipy
try:
    import scipy  # noqa
except ImportError:
    moma = None
else:
    from cobra.flux_analysis import moma


def single_reaction_deletion(cobra_model, reaction_list=None, solver=None,
                             method="fba", **solver_args):
    """Sequentially knocks out each reaction from a given reaction list.

    Parameters
    ----------
    cobra_model : a cobra model
        The model from which to delete the reactions. The model will not be
        modified.
    reaction_list : iterable
        List of reaction IDs or cobra.Reaction. If None (default) will use all
        reactions in the model.
    method : str, optional
        The method used to obtain fluxes. Must be one of "fba" or "moma".
    solver : str, optional
        Name of the solver to be used.
    method : str, optional
        The method used to obtain fluxes. Must be one of "fba" or "moma".
    solver_args : optional
        Additional arguments for the solver. Ignored for optlang solver, please
        use `model.solver.configuration` instead.

    Returns
    -------
    pandas.DataFrame
        Data frame with two column and reaction id as index:
        - flux: the value of the objective after the knockout
        - status: the solution's status, (for instance "optimal" for each
          knockout)
    """
    if reaction_list is None:
        reaction_list = cobra_model.reactions
    else:
        reaction_list = [cobra_model.reactions.get_by_id(i)
                         if isinstance(i, string_types) else i
                         for i in reaction_list]
    if method == "fba":
        result = single_reaction_deletion_fba(cobra_model, reaction_list,
                                              solver=solver, **solver_args)
    elif method == "moma":
        result = single_reaction_deletion_moma(cobra_model, reaction_list,
                                               solver=solver, **solver_args)
    else:
        raise ValueError("Unknown deletion method '%s'" % method)
    return pandas.DataFrame({'flux': result[0], 'status': result[1]})


def single_reaction_deletion_fba(cobra_model, reaction_list, solver=None,
                                 **solver_args):
    """Sequentially knocks out each reaction in a model using FBA.

    Not supposed to be called directly use
    `single_reactions_deletion(..., method="fba")` instead.

    Parameters
    ----------
    reaction_list : iterable
        List of reaction Ids or cobra.Reaction.
    solver: str, optional
        The name of the solver to be used.

    Returns
    -------
    tuple of dicts
        A tuple ({reaction_id: growth_rate}, {reaction_id: status})
    """
    legacy = False
    if solver is None:
        solver = cobra_model.solver
    elif "optlang-" in solver:
        solver = solvers.interface_to_str(solver)
        solver = solvers.solvers[solver]
    else:
        legacy = True
        solver = legacy_solvers.solver_dict[solver]
        lp = solver.create_problem(cobra_model)

    growth_rate_dict = {}
    status_dict = {}

    if not legacy:
        with cobra_model as m:
            m.solver = solver
            for reaction in reaction_list:
                with m:
                    reaction.bounds = (0.0, 0.0)
                    m.solver.optimize()
                    status = m.solver.status
                    status_dict[reaction.id] = status
                    growth_rate_dict[reaction.id] = m.solver.objective.value \
                        if status == OPTIMAL else 0.
    else:
        # This entire block can be removed once the legacy solvers are
        # deprecated
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
            solver.change_variable_bounds(lp, index, old_bounds[0],
                                          old_bounds[1])
    return growth_rate_dict, status_dict


def single_reaction_deletion_moma(cobra_model, reaction_list, solver=None,
                                  **solver_args):
    """Sequentially knocks out each reaction in a model using MOMA.

    Not supposed to be called directly use
    `single_reactions_deletion(..., method="moma")` instead.

    Parameters
    ----------
    reaction_list : iterable
        List of reaction IDs or cobra.Reaction.
    solver: str, optional
        The name of the solver to be used.

    Returns
    -------
    tuple of dicts
        A tuple ({reaction_id: growth_rate}, {reaction_id: status})
    """
    # The same function can not be used because MOMA can not re-use the
    # same LP object. Problem re-use leads to incorrect solutions.
    # This is *not* true for optlang solvers!
    if moma is None:
        raise RuntimeError("scipy required for moma")

    legacy = False
    if solver is None:
        solver = cobra_model.solver
    elif "optlang-" in solver:
        solver = solvers.interface_to_str(solver)
        solver = solvers.solvers[solver]
    else:
        legacy = True
        solver = legacy_solvers.solver_dict[solver]
        moma_model, moma_objective = moma.create_euclidian_moma_model(
            cobra_model)

    growth_rate_dict = {}
    status_dict = {}

    if not legacy:
        with cobra_model as m:
            m.solver = solver
            moma.add_moma(m)
            for reaction in reaction_list:
                with m:
                    reaction.bounds = (0.0, 0.0)
                    m.solver.optimize()
                    status = m.solver.status
                    status_dict[reaction.id] = status
                    if status == OPTIMAL:
                        growth = m.variables.moma_old_objective.primal
                    else:
                        growth = 0.0
                    growth_rate_dict[reaction.id] = growth
    else:
        for reaction in reaction_list:
            index = cobra_model.reactions.index(reaction)
            solution = moma.moma_knockout(moma_model, moma_objective, (index,),
                                          solver=solver, **solver_args)
            status_dict[reaction.id] = solution.status
            growth_rate_dict[reaction.id] = solution.f
    return growth_rate_dict, status_dict


def single_gene_deletion(cobra_model, gene_list=None, solver=None,
                         method="fba", **solver_args):
    """Sequentially knocks out each gene from a given gene list.

    Parameters
    ----------
    cobra_model : a cobra model
        The model from which to delete the genes. The model will not be
        modified.
    gene_list : iterable
        List of gene IDs or cobra.Gene. If None (default) will use all genes in
        the model.
    method : str, optional
        The method used to obtain fluxes. Must be one of "fba" or "moma".
    solver : str, optional
        Name of the solver to be used.
    solver_args : optional
        Additional arguments for the solver. Ignored for optlang solver, please
        use `model.solver.configuration` instead.

    Returns
    -------
    pandas.DataFrame
        Data frame with two column and reaction id as index:
        - flux: the value of the objective after the knockout
        - status: the solution's status, (for instance "optimal" for each
          knockout)
    """
    if gene_list is None:
        gene_list = cobra_model.genes
    else:
        gene_list = [cobra_model.genes.get_by_id(i)
                     if isinstance(i, string_types) else i for i in gene_list]

    if method == "fba":
        result = single_gene_deletion_fba(cobra_model, gene_list,
                                          solver=solver, **solver_args)
    elif method == "moma":
        result = single_gene_deletion_moma(cobra_model, gene_list,
                                           solver=solver, **solver_args)
    else:
        raise ValueError("Unknown deletion method '%s'" % method)
    return pandas.DataFrame({'flux': result[0], 'status': result[1]})


def single_gene_deletion_fba(cobra_model, gene_list, solver=None,
                             **solver_args):
    """Sequentially knocks out each gene in a model using FBA.

    Not supposed to be called directly use
    `single_reactions_deletion(..., method="fba")` instead.

    Parameters
    ----------
    gene_list : iterable
        List of gene IDs or cobra.Reaction.
    solver: str, optional
        The name of the solver to be used.

    Returns
    -------
    tuple of dicts
        A tuple ({reaction_id: growth_rate}, {reaction_id: status})
    """
    legacy = False
    if solver is None:
        solver = cobra_model.solver
    elif "optlang-" in solver:
        solver = solvers.interface_to_str(solver)
        solver = solvers.solvers[solver]
    else:
        legacy = True
        solver = legacy_solvers.solver_dict[solver]
        lp = solver.create_problem(cobra_model)

    growth_rate_dict = {}
    status_dict = {}

    if not legacy:
        with cobra_model as m:
            m.solver = solver
            for gene in gene_list:
                ko = find_gene_knockout_reactions(cobra_model, [gene])
                with m:
                    for reaction in ko:
                        reaction.bounds = (0.0, 0.0)
                    m.solver.optimize()
                    status = m.solver.status
                    status_dict[gene.id] = status
                    growth_rate_dict[gene.id] = m.solver.objective.value if \
                        status == OPTIMAL else 0.
    else:
        for gene in gene_list:
            old_bounds = {}
            for reaction in find_gene_knockout_reactions(cobra_model, [gene]):
                index = cobra_model.reactions.index(reaction)
                old_bounds[index] = reaction.bounds
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
    return growth_rate_dict, status_dict


def single_gene_deletion_moma(cobra_model, gene_list, solver=None,
                              **solver_args):
    """Sequentially knocks out each gene in a model using MOMA.

    Not supposed to be called directly use
    `single_reactions_deletion(..., method="moma")` instead.

    Parameters
    ----------
    gene_list : iterable
        List of gene IDs or cobra.Reaction.
    solver: str, optional
        The name of the solver to be used.

    Returns
    -------
    tuple of dicts
        A tuple ({reaction_id: growth_rate}, {reaction_id: status})
    """
    if moma is None:
        raise RuntimeError("scipy required for moma")

    legacy = False
    if solver is None:
        solver = cobra_model.solver
    elif "optlang-" in solver:
        solver = solvers.interface_to_str(solver)
        solver = solvers.solvers[solver]
    else:
        legacy = True
        solver = legacy_solvers.solver_dict[solver]
        moma_model, moma_objective = moma.create_euclidian_moma_model(
            cobra_model)

    growth_rate_dict = {}
    status_dict = {}

    if not legacy:
        with cobra_model as m:
            m.solver = solver
            moma.add_moma(m)
            for gene in gene_list:
                ko = find_gene_knockout_reactions(cobra_model, [gene])
                with m:
                    for reaction in ko:
                        reaction.bounds = (0.0, 0.0)
                    m.solver.optimize()
                    status = m.solver.status
                    status_dict[gene.id] = status
                    if status == "optimal":
                        growth = m.variables.moma_old_objective.primal
                    else:
                        growth = 0.0
                    growth_rate_dict[gene.id] = growth
    else:
        for gene in gene_list:
            delete_model_genes(moma_model, [gene.id])
            solution = moma.solve_moma_model(moma_model, moma_objective,
                                             solver=solver, **solver_args)
            status_dict[gene.id] = solution.status
            growth_rate_dict[gene.id] = solution.f
            undelete_model_genes(moma_model)
    return growth_rate_dict, status_dict
