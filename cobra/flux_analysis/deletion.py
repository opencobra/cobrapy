# -*- coding: utf-8 -*-

import multiprocessing
import logging
import optlang
from warnings import warn
from itertools import product
from functools import partial
from builtins import (map, dict)

import pandas as pd

from cobra.manipulation.delete import find_gene_knockout_reactions
import cobra.util.solver as sutil
from cobra.flux_analysis.moma import add_moma

LOGGER = logging.getLogger(__name__)


def _reactions_knockouts_with_restore(model, reactions):
    with model:
        for reaction in reactions:
            reaction.knock_out()
        growth = _get_growth(model)
    return [r.id for r in reactions], growth, model.solver.status


def _get_growth(model):
    try:
        if 'moma_old_objective' in model.solver.variables:
            model.slim_optimize()
            growth = model.solver.variables.moma_old_objective.primal
        else:
            growth = model.slim_optimize()
    except optlang.exceptions.SolverError:
        growth = float('nan')
    return growth


def _reaction_deletion(model, ids):
    return _reactions_knockouts_with_restore(
        model,
        [model.reactions.get_by_id(r_id) for r_id in ids]
    )


def _gene_deletion(model, ids):
    all_reactions = []
    for g_id in ids:
        all_reactions.extend(
            find_gene_knockout_reactions(
                model, (model.genes.get_by_id(g_id),)
            )
        )
    _, growth, status = _reactions_knockouts_with_restore(model, all_reactions)
    return (ids, growth, status)


def _reaction_deletion_worker(ids):
    global _model
    return _reaction_deletion(_model, ids)


def _gene_deletion_worker(ids):
    global _model
    return _gene_deletion(_model, ids)


def _init_worker(model):
    global _model
    _model = model


def _multi_deletion(model, entity, element_lists, method="fba",
                    processes=None):
    """
    Provide a common interface for single or multiple knockouts.

    Parameters
    ----------
    model : cobra.Model
        The metabolic model to perform deletions in.

    entity : 'gene' or 'reaction'
        The entity to knockout (``cobra.Gene`` or ``cobra.Reaction``).

    element_lists : list
        List of iterables ``cobra.Reaction``s or ``cobra.Gene``s (or their IDs)
        to be deleted.

    method: {"fba", "moma", "linear moma"}, optional
        Method used to predict the growth rate.

    processes : int, optional
        The number of parallel processes to run. Can speed up the computations
        if the number of knockouts to perform is large. If not passed,
        will be set to the number of CPUs found.

    Returns
    -------
    pandas.DataFrame
        A representation of all combinations of entity deletions. The
        columns are 'growth' and 'status', where

        index : frozenset([str])
            The gene or reaction identifiers that were knocked out.
        growth : float
            The growth rate of the adjusted model.
        status : str
            The solution's status.
    """
    solver = sutil.interface_to_str(model.problem.__name__)
    if "moma" in method and solver not in sutil.qp_solvers:
        raise RuntimeError(
            "Cannot use MOMA since '{}' is not QP-capable."
            "Please choose a different solver or use FBA only.".format(solver))

    if processes is None:
        try:
            processes = multiprocessing.cpu_count()
        except NotImplementedError:
            warn("Number of cores could not be detected - assuming 1.")
            processes = 1

    with model:
        if "moma" in method:
            add_moma(model, linear="linear" in method)

        args = set([frozenset(comb) for comb in product(*element_lists)])
        processes = min(processes, len(args))

        def extract_knockout_results(result_iter):
            result = pd.DataFrame([
                (frozenset(ids), growth, status)
                for (ids, growth, status) in result_iter
            ], columns=['ids', 'growth', 'status'])
            result.set_index('ids', inplace=True)
            return result

        if processes > 1:
            worker = dict(gene=_gene_deletion_worker,
                          reaction=_reaction_deletion_worker)[entity]
            chunk_size = len(args) // processes
            pool = multiprocessing.Pool(
                processes, initializer=_init_worker, initargs=(model,)
            )
            results = extract_knockout_results(pool.imap_unordered(
                worker,
                args,
                chunksize=chunk_size
            ))
            pool.close()
            pool.join()
        else:
            worker = dict(gene=_gene_deletion,
                          reaction=_reaction_deletion)[entity]
            results = extract_knockout_results(map(
                partial(worker, model), args))
        return results


def _entities_ids(entities):
    try:
        return [e.id for e in entities]
    except AttributeError:
        return list(entities)


def _element_lists(entities, *ids):
    lists = list(ids)
    if lists[0] is None:
        lists[0] = entities
    result = [_entities_ids(lists[0])]
    for l in lists[1:]:
        if l is None:
            result.append(result[-1])
        else:
            result.append(_entities_ids(l))
    return result


def single_reaction_deletion(model, reaction_list=None, method="fba",
                             processes=None):
    """
    Knock out each reaction from a given list.

    Parameters
    ----------
    model : cobra.Model
        The metabolic model to perform deletions in.

    reaction_list : iterable
        ``cobra.Reaction``s to be deleted. If not passed,
        all the reactions from the model are used.

    method: {"fba", "moma", "linear moma"}, optional
        Method used to predict the growth rate.

    processes : int, optional
        The number of parallel processes to run. Can speed up the computations
        if the number of knockouts to perform is large. If not passed,
        will be set to the number of CPUs found.

    Returns
    -------
    pandas.DataFrame
        A representation of all single reaction deletions. The columns are
        'growth' and 'status', where

        index : frozenset([str])
            The reaction identifier that was knocked out.
        growth : float
            The growth rate of the adjusted model.
        status : str
            The solution's status.

    """
    return _multi_deletion(
        model, 'reaction',
        element_lists=_element_lists(model.reactions, reaction_list),
        method=method, processes=processes)


def single_gene_deletion(model, gene_list=None, method="fba", processes=None):
    """
    Knock out each gene from a given list.

    Parameters
    ----------
    model : cobra.Model
        The metabolic model to perform deletions in.

    gene_list : iterable
        ``cobra.Gene``s to be deleted. If not passed,
        all the genes from the model are used.

    method: {"fba", "moma", "linear moma"}, optional
        Method used to predict the growth rate.

    processes : int, optional
        The number of parallel processes to run. Can speed up the computations
        if the number of knockouts to perform is large. If not passed,
        will be set to the number of CPUs found.

    Returns
    -------
    pandas.DataFrame
        A representation of all single gene deletions. The columns are
        'growth' and 'status', where

        index : frozenset([str])
            The gene identifier that was knocked out.
        growth : float
            The growth rate of the adjusted model.
        status : str
            The solution's status.

    """
    return _multi_deletion(
        model, 'gene', element_lists=_element_lists(model.genes, gene_list),
        method=method, processes=processes)


def double_reaction_deletion(model, reaction_list1=None, reaction_list2=None,
                             method="fba", processes=None):
    """
    Knock out each reaction pair from the combinations of two given lists.

    We say 'pair' here but the order order does not matter.

    Parameters
    ----------
    model : cobra.Model
        The metabolic model to perform deletions in.

    reaction_list1 : iterable, optional
        First iterable of ``cobra.Reaction``s to be deleted. If not passed,
        all the reactions from the model are used.

    reaction_list2 : iterable, optional
        Second iterable of ``cobra.Reaction``s to be deleted. If not passed,
        all the reactions from the model are used.

    method: {"fba", "moma", "linear moma"}, optional
        Method used to predict the growth rate.

    processes : int, optional
        The number of parallel processes to run. Can speed up the computations
        if the number of knockouts to perform is large. If not passed,
        will be set to the number of CPUs found.

    Returns
    -------
    pandas.DataFrame
        A representation of all combinations of reaction deletions. The
        columns are 'growth' and 'status', where

        index : frozenset([str])
            The reaction identifiers that were knocked out.
        growth : float
            The growth rate of the adjusted model.
        status : str
            The solution's status.

    """

    reaction_list1, reaction_list2 = _element_lists(model.reactions,
                                                    reaction_list1,
                                                    reaction_list2)
    return _multi_deletion(
        model, 'reaction', element_lists=[reaction_list1, reaction_list2],
        method=method, processes=processes)


def double_gene_deletion(model, gene_list1=None, gene_list2=None,
                         method="fba", processes=None):
    """
    Knock out each gene pair from the combination of two given lists.

    We say 'pair' here but the order order does not matter.

    Parameters
    ----------
    model : cobra.Model
        The metabolic model to perform deletions in.

    gene_list1 : iterable, optional
        First iterable of ``cobra.Gene``s to be deleted. If not passed,
        all the genes from the model are used.

    gene_list2 : iterable, optional
        Second iterable of ``cobra.Gene``s to be deleted. If not passed,
        all the genes from the model are used.

    method: {"fba", "moma", "linear moma"}, optional
        Method used to predict the growth rate.

    processes : int, optional
        The number of parallel processes to run. Can speed up the computations
        if the number of knockouts to perform is large. If not passed,
        will be set to the number of CPUs found.

    Returns
    -------
    pandas.DataFrame
        A representation of all combinations of gene deletions. The
        columns are 'growth' and 'status', where

        index : frozenset([str])
            The gene identifiers that were knocked out.
        growth : float
            The growth rate of the adjusted model.
        status : str
            The solution's status.

    """

    gene_list1, gene_list2 = _element_lists(model.genes, gene_list1,
                                            gene_list2)
    return _multi_deletion(
        model, 'gene', element_lists=[gene_list1, gene_list2],
        method=method, processes=processes)
