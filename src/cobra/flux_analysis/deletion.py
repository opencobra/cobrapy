# -*- coding: utf-8 -*-

import logging
import multiprocessing
from builtins import dict, map
from functools import partial
from itertools import product
from typing import List, Set, Union

import pandas as pd
from optlang.exceptions import SolverError

from cobra.core import Configuration, Gene, Reaction
from cobra.flux_analysis.moma import add_moma
from cobra.flux_analysis.room import add_room
from cobra.manipulation.delete import find_gene_knockout_reactions
from cobra.util import solver as sutil


LOGGER = logging.getLogger(__name__)
CONFIGURATION = Configuration()


def _reactions_knockouts_with_restore(model, reactions):
    with model:
        for reaction in reactions:
            reaction.knock_out()
        growth = _get_growth(model)
    return [r.id for r in reactions], growth, model.solver.status


def _get_growth(model):
    try:
        if "moma_old_objective" in model.solver.variables:
            model.slim_optimize()
            growth = model.solver.variables.moma_old_objective.primal
        else:
            growth = model.slim_optimize()
    except SolverError:
        growth = float("nan")
    return growth


def _reaction_deletion(model, ids):
    return _reactions_knockouts_with_restore(
        model, [model.reactions.get_by_id(r_id) for r_id in ids]
    )


def _gene_deletion(model, ids):
    all_reactions = []
    for g_id in ids:
        all_reactions.extend(
            find_gene_knockout_reactions(model, (model.genes.get_by_id(g_id),))
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


def _multi_deletion(
    model, entity, element_lists, method="fba", solution=None, processes=None, **kwargs
):
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
    method: {"fba", "moma", "linear moma", "room", "linear room"}, optional
        Method used to predict the growth rate.
    solution : cobra.Solution, optional
        A previous solution to use as a reference for (linear) MOMA or ROOM.
    processes : int, optional
        The number of parallel processes to run. Can speed up the computations
        if the number of knockouts to perform is large. If not passed,
        will be set to the number of CPUs found.
    kwargs :
        Passed on to underlying simulation functions.

    Returns
    -------
    pandas.DataFrame
        A representation of all combinations of entity deletions. The
        columns are 'growth' and 'status', where

        index : tuple(str)
            The gene or reaction identifiers that were knocked out.
        growth : float
            The growth rate of the adjusted model.
        status : str
            The solution's status.
    """
    solver = sutil.interface_to_str(model.problem.__name__)
    if method == "moma" and solver not in sutil.qp_solvers:
        raise RuntimeError(
            "Cannot use MOMA since '{}' is not QP-capable."
            "Please choose a different solver or use FBA only.".format(solver)
        )

    if processes is None:
        processes = CONFIGURATION.processes

    with model:
        if "moma" in method:
            add_moma(model, solution=solution, linear="linear" in method)
        elif "room" in method:
            add_room(model, solution=solution, linear="linear" in method, **kwargs)

        args = set([frozenset(comb) for comb in product(*element_lists)])
        processes = min(processes, len(args))

        def extract_knockout_results(result_iter):
            result = pd.DataFrame(
                [(set(ids), growth, status,) for (ids, growth, status) in result_iter],
                columns=["ids", "growth", "status"],
            )
            return result

        if processes > 1:
            worker = dict(
                gene=_gene_deletion_worker, reaction=_reaction_deletion_worker
            )[entity]
            chunk_size = len(args) // processes
            pool = multiprocessing.Pool(
                processes, initializer=_init_worker, initargs=(model,)
            )
            results = extract_knockout_results(
                pool.imap_unordered(worker, args, chunksize=chunk_size)
            )
            pool.close()
            pool.join()
        else:
            worker = dict(gene=_gene_deletion, reaction=_reaction_deletion)[entity]
            results = extract_knockout_results(map(partial(worker, model), args))
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


def single_reaction_deletion(
    model, reaction_list=None, method="fba", solution=None, processes=None, **kwargs
):
    """
    Knock out each reaction from a given list.

    Parameters
    ----------
    model : cobra.Model
        The metabolic model to perform deletions in.
    reaction_list : iterable, optional
        ``cobra.Reaction``s to be deleted. If not passed,
        all the reactions from the model are used.
    method: {"fba", "moma", "linear moma", "room", "linear room"}, optional
        Method used to predict the growth rate.
    solution : cobra.Solution, optional
        A previous solution to use as a reference for (linear) MOMA or ROOM.
    processes : int, optional
        The number of parallel processes to run. Can speed up the computations
        if the number of knockouts to perform is large. If not passed,
        will be set to the number of CPUs found.
    kwargs :
        Keyword arguments are passed on to underlying simulation functions
        such as ``add_room``.

    Returns
    -------
    pandas.DataFrame
        A representation of all single reaction deletions. The columns are
        'growth' and 'status', where

        index : tuple(str)
            The reaction identifier that was knocked out.
        growth : float
            The growth rate of the adjusted model.
        status : str
            The solution's status.

    """
    return _multi_deletion(
        model,
        "reaction",
        element_lists=_element_lists(model.reactions, reaction_list),
        method=method,
        solution=solution,
        processes=processes,
        **kwargs
    )


def single_gene_deletion(
    model, gene_list=None, method="fba", solution=None, processes=None, **kwargs
):
    """
    Knock out each gene from a given list.

    Parameters
    ----------
    model : cobra.Model
        The metabolic model to perform deletions in.
    gene_list : iterable
        ``cobra.Gene``s to be deleted. If not passed,
        all the genes from the model are used.
    method: {"fba", "moma", "linear moma", "room", "linear room"}, optional
        Method used to predict the growth rate.
    solution : cobra.Solution, optional
        A previous solution to use as a reference for (linear) MOMA or ROOM.
    processes : int, optional
        The number of parallel processes to run. Can speed up the computations
        if the number of knockouts to perform is large. If not passed,
        will be set to the number of CPUs found.
    kwargs :
        Keyword arguments are passed on to underlying simulation functions
        such as ``add_room``.

    Returns
    -------
    pandas.DataFrame
        A representation of all single gene deletions. The columns are
        'growth' and 'status', where

        index : tuple(str)
            The gene identifier that was knocked out.
        growth : float
            The growth rate of the adjusted model.
        status : str
            The solution's status.

    """
    return _multi_deletion(
        model,
        "gene",
        element_lists=_element_lists(model.genes, gene_list),
        method=method,
        solution=solution,
        processes=processes,
        **kwargs
    )


def double_reaction_deletion(
    model,
    reaction_list1=None,
    reaction_list2=None,
    method="fba",
    solution=None,
    processes=None,
    **kwargs
):
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
    method: {"fba", "moma", "linear moma", "room", "linear room"}, optional
        Method used to predict the growth rate.
    solution : cobra.Solution, optional
        A previous solution to use as a reference for (linear) MOMA or ROOM.
    processes : int, optional
        The number of parallel processes to run. Can speed up the computations
        if the number of knockouts to perform is large. If not passed,
        will be set to the number of CPUs found.
    kwargs :
        Keyword arguments are passed on to underlying simulation functions
        such as ``add_room``.

    Returns
    -------
    pandas.DataFrame
        A representation of all combinations of reaction deletions. The
        columns are 'growth' and 'status', where

        index : tuple(str)
            The reaction identifiers that were knocked out.
        growth : float
            The growth rate of the adjusted model.
        status : str
            The solution's status.

    """

    reaction_list1, reaction_list2 = _element_lists(
        model.reactions, reaction_list1, reaction_list2
    )
    return _multi_deletion(
        model,
        "reaction",
        element_lists=[reaction_list1, reaction_list2],
        method=method,
        solution=solution,
        processes=processes,
        **kwargs
    )


def double_gene_deletion(
    model,
    gene_list1=None,
    gene_list2=None,
    method="fba",
    solution=None,
    processes=None,
    **kwargs
):
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
    method: {"fba", "moma", "linear moma", "room", "linear room"}, optional
        Method used to predict the growth rate.
    solution : cobra.Solution, optional
        A previous solution to use as a reference for (linear) MOMA or ROOM.
    processes : int, optional
        The number of parallel processes to run. Can speed up the computations
        if the number of knockouts to perform is large. If not passed,
        will be set to the number of CPUs found.
    kwargs :
        Keyword arguments are passed on to underlying simulation functions
        such as ``add_room``.

    Returns
    -------
    pandas.DataFrame
        A representation of all combinations of gene deletions. The
        columns are 'growth' and 'status', where

        index : tuple(str)
            The gene identifiers that were knocked out.
        growth : float
            The growth rate of the adjusted model.
        status : str
            The solution's status.

    """

    gene_list1, gene_list2 = _element_lists(model.genes, gene_list1, gene_list2)
    return _multi_deletion(
        model,
        "gene",
        element_lists=[gene_list1, gene_list2],
        method=method,
        solution=solution,
        processes=processes,
        **kwargs
    )


@pd.api.extensions.register_dataframe_accessor("knockout")
class KnockoutAccessor:
    """Access unique combinations of reactions in deletion results.

    This allows acces in the form of `results.knockout[rxn1]` or
    `results.knockout["rxn1_id"]`. Each individual entry will return a deletion
    so `results.knockout[rxn1, rxn2]` will return two deletions (for individual
    knockouts of rxn1 and rxn2 respectively). Multi-deletions can be accessed by passing
    in sets like `results.knockou[{rxn1, rxn2}]` which denotes the double deletion of
    both reactions. Thus, the following are allowed index elements:

    - single reactions or genes (depending on whether it is a gene or reaction deletion)
    - single reaction IDs or gene IDs
    - lists of single single reaction IDs or gene IDs (will return one row for each
      element in the list)
    - sets of reactions or genes (for multi-deletions)
    - sets of reactions IDs or gene IDs
    - list of sets of objects or IDs (to get several multi-deletions)
    """

    def __init__(self, pandas_obj: pd.DataFrame) -> None:
        """Set up the accessor.

        Parameters:
        -----------
        pandas_obj : pd.DataFrame or pd.Series
            A result from one of the deletion methods.
        """
        self._validate(pandas_obj)
        self._result = pandas_obj

    @staticmethod
    def _validate(obj: pd.DataFrame) -> None:
        # verify it is a deletion results
        if any(name not in obj.columns for name in ["ids", "growth", "status"]):
            raise AttributeError("Must be DataFrame returned by a deletion method.")

    def __getitem__(
        self,
        args: Union[
            Gene,
            List[Gene],
            Set[Gene],
            List[Set[Gene]],
            Reaction,
            List[Reaction],
            Set[Reaction],
            List[Set[Reaction]],
            str,
            List[str],
            Set[str],
            List[Set[str]],
        ],
    ) -> pd.DataFrame:
        """Return the deletion result for a particular set of knocked entities.

        Parameters:
        -----------
        args : cobra.Reactions, cobra.Gene, str, set, or list
            The deletions to be returned. Accepts:
            - single reactions or genes
            - single reaction IDs or gene IDs
            - lists of single single reaction IDs or gene IDs
            - sets of reactions or genes
            - sets of reactions IDs or gene IDs
            - list of sets of objects or IDs
            See the docs for usage examples.

        Returns:
        --------
        pd.DataFrame
            The deletion result where the chosen entities have been deleted. Each row
            denotes a deletion.
        """
        if not any(isinstance(args, t) for t in [tuple, list]):
            args = [args]

        if any(isinstance(args[0], t) for t in [Reaction, Gene, str]):
            try:
                args = [{obj.id} for obj in args]
            except AttributeError:
                # are already strings
                args = [{obj} for obj in args]
        elif isinstance(args[0], set):
            try:
                args = [set(elem.id for elem in obj) for obj in args]
            except AttributeError:
                args = [set(obj) for obj in args]
        else:
            raise ValueError(
                "Allowed indices are single Reactions or Genes, "
                "lists of Reactions of Genes, or lists of sets "
                "of Reactions or Genes."
            )
        found = [x in args for x in self._result.ids]
        return self._result[found]
