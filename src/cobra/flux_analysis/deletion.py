"""Provide functions for reaction and gene deletions."""

from functools import partial
from itertools import product
from typing import TYPE_CHECKING, List, Optional, Set, Tuple, Union

import pandas as pd
from optlang.exceptions import SolverError

from ..core import Configuration, Gene, Model, Reaction
from ..util import ProcessPool
from ..util import solver as sutil
from .moma import add_moma
from .room import add_room


if TYPE_CHECKING:
    from cobra import Solution


configuration = Configuration()


def _get_growth(model: Model) -> Tuple[float, str]:
    """Return the growth from the `model`.

    Parameters
    ----------
    model : cobra.Model
        The model to obtain growth for.

    Returns
    -------
    float
        The obtained growth value. Returns nan if there is some error while
        optimizing.

    """
    try:
        if "moma_old_objective" in model.solver.variables:
            model.slim_optimize()
            growth = model.solver.variables.moma_old_objective.primal
        else:
            growth = model.slim_optimize()
    except SolverError:
        growth = float("nan")
    return growth, model.solver.status


def _reaction_deletion(
    model: Model, reaction_ids: List[str]
) -> Tuple[List[str], float, str]:
    """Perform reaction deletion.

    Parameters
    ----------
    model : cobra.Model
        The model to perform reaction deletion on.
    ids : list of str
        The reaction IDs to knock-out.

    Returns
    -------
    tuple of (list of str, float, str)
        A tuple containing reaction IDs knocked out, growth of the model and
        the solver status.

    """
    with model:
        for rxn_id in reaction_ids:
            model.reactions.get_by_id(rxn_id).knock_out()
        growth, status = _get_growth(model)
    return reaction_ids, growth, status


def _reaction_deletion_worker(ids: List[str]) -> Tuple[List[str], float, str]:
    """Perform reaction deletions on worker process.

    Parameters
    ----------
    ids : list of str
        The reaction IDs to knock-out.

    Returns
    -------
    tuple of (list of str, float, str)
        A tuple containing reaction IDs knocked out, growth of the model and
        the solver status.

    """
    global _model

    return _reaction_deletion(_model, ids)


def _gene_deletion(model: Model, gene_ids: List[str]) -> Tuple[List[str], float, str]:
    """Perform gene deletions.

    Parameters
    ----------
    model : cobra.Model
        The model to perform gene deletion on.
    ids : list of str
        The gene IDs to knock-out.

    Returns
    -------
    tuple of (list of str, float, str)
        A tuple containing gene IDs knocked out, growth of the model and
        the solver status.

    """
    with model:
        for gene_id in gene_ids:
            model.genes.get_by_id(gene_id).knock_out()
        growth, status = _get_growth(model)
    return gene_ids, growth, status


def _gene_deletion_worker(ids: List[str]) -> Tuple[List[str], float, str]:
    """Perform gene deletions on worker process.

    Parameters
    ----------
    ids : list of str
        The gene IDs to knock-out.

    Returns
    -------
    tuple of (list of str, float, str)
        A tuple containing gene IDs knocked out, growth of the model and
        the solver status.

    """
    global _model

    return _gene_deletion(_model, ids)


def _init_worker(model: Model) -> None:
    """Initialize worker process."""
    global _model

    _model = model


def _multi_deletion(
    model: Model,
    entity: str,
    element_lists: List[Union[Gene, Reaction]],
    method: str = "fba",
    solution: Optional["Solution"] = None,
    processes: Optional[int] = None,
    **kwargs,
) -> pd.DataFrame:
    """Provide a common interface for single or multiple knockouts.

    Parameters
    ----------
    model : cobra.Model
        The metabolic model to perform deletions in.
    entity : {"gene", "reaction"}
        The entity to knockout.
    element_lists : list of cobra.Gene or cobra.Reaction
        List of cobra.Gene or cobra.Reaction to be deleted.
    method : {"fba", "moma", "linear moma", "room", "linear room"}, optional
        Method used to predict the growth rate (default "fba").
    solution : cobra.Solution, optional
        A previous solution to use as a reference for (linear) MOMA or ROOM
        (default None).
    processes : int, optional
        The number of parallel processes to run. Can speed up the computations
        if the number of knockouts to perform is large. If not passed,
        will be set to `configuration.processes` (default None).
    **kwargs :
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
            f"Cannot use MOMA since '{solver}' is not QP-capable. "
            "Please choose a different solver or use FBA only."
        )

    if processes is None:
        processes = configuration.processes

    with model:
        if "moma" in method:
            add_moma(model, solution=solution, linear="linear" in method)
        elif "room" in method:
            add_room(model, solution=solution, linear="linear" in method, **kwargs)

        args = {frozenset(comb) for comb in product(*element_lists)}
        processes = min(processes, len(args))

        def extract_knockout_results(result_iter):
            result = pd.DataFrame(
                [
                    (
                        set(ids),
                        growth,
                        status,
                    )
                    for (ids, growth, status) in result_iter
                ],
                columns=["ids", "growth", "status"],
            )
            return result

        if processes > 1:
            worker = {
                "gene": _gene_deletion_worker,
                "reaction": _reaction_deletion_worker,
            }[entity]
            chunk_size = len(args) // processes

            with ProcessPool(
                processes, initializer=_init_worker, initargs=(model,)
            ) as pool:
                results = extract_knockout_results(
                    pool.imap_unordered(worker, args, chunksize=chunk_size)
                )
        else:
            worker = {"gene": _gene_deletion, "reaction": _reaction_deletion}[entity]
            results = extract_knockout_results(map(partial(worker, model), args))
        return results


def _entities_ids(entities: List[Union[str, Gene, Reaction]]) -> List[str]:
    """Return the IDs of the `entities`.

    Parameters
    ----------
    entities : list of str or cobra.Gene or cobra.Reaction
        The list of entities whose IDs need to be returned.

    Returns
    -------
    list of str
        The IDs of the `entities`.

    """
    try:
        return [e.id for e in entities]
    except AttributeError:
        return list(entities)


def _element_lists(
    entities: List[Union[str, Gene, Reaction]], *ids: List[str]
) -> List[str]:
    """Return the elements.

    Parameters
    ----------
    entities : list of str or cobra.Gene or cobra.Reaction
        The list of entities.
    *ids : list of str
        The list of IDs.

    Returns
    -------
    list of str
        The list of IDs.

    """
    lists = list(ids)
    if lists[0] is None:
        lists[0] = entities
    result = [_entities_ids(lists[0])]
    for _list in lists[1:]:
        if _list is None:
            result.append(result[-1])
        else:
            result.append(_entities_ids(_list))
    return result


def single_reaction_deletion(
    model: Model,
    reaction_list: Optional[List[Union[Reaction, str]]] = None,
    method: str = "fba",
    solution: Optional["Solution"] = None,
    processes: Optional[int] = None,
    **kwargs,
) -> pd.DataFrame:
    """Knock out each reaction from `reaction_list`.

    Parameters
    ----------
    model : cobra.Model
        The metabolic model to perform deletions in.
    reaction_list : list of cobra.Reaction or str, optional
        The reactions be knocked out. If not passed, all the reactions from
        the model are used (default None).
    method: {"fba", "moma", "linear moma", "room", "linear room"}, optional
        Method used to predict the growth rate (default "fba").
    solution : cobra.Solution, optional
        A previous solution to use as a reference for (linear) MOMA or ROOM
        (default None).
    processes : int, optional
        The number of parallel processes to run. Can speed up the computations
        if the number of knockouts to perform is large. If not passed,
        will be set to `configuration.processes` (default None).
    **kwargs :
        Keyword arguments are passed on to underlying simulation functions
        such as `add_room`.

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
        **kwargs,
    )


def single_gene_deletion(
    model: Model,
    gene_list: Optional[List[Union[Gene, str]]] = None,
    method: str = "fba",
    solution: Optional["Solution"] = None,
    processes: Optional[int] = None,
    **kwargs,
) -> pd.DataFrame:
    """Knock out each gene from `gene_list`.

    Parameters
    ----------
    model : cobra.Model
        The metabolic model to perform deletions in.
    gene_list : list of cobra.Gene or str, optional
        The gene objects to be deleted. If not passed, all the genes from the
        model are used (default None).
    method : {"fba", "moma", "linear moma", "room", "linear room"}, optional
        Method used to predict the growth rate (default "fba").
    solution : cobra.Solution, optional
        A previous solution to use as a reference for (linear) MOMA or ROOM
        (default None).
    processes : int, optional
        The number of parallel processes to run. Can speed up the computations
        if the number of knockouts to perform is large. If not passed,
        will be set to `configuration.processes` (default None).
    **kwargs :
        Keyword arguments are passed on to underlying simulation functions
        such as `add_room`.

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
        **kwargs,
    )


def double_reaction_deletion(
    model: Model,
    reaction_list1: Optional[List[Union[Reaction, str]]] = None,
    reaction_list2: Optional[List[Union[Reaction, str]]] = None,
    method: str = "fba",
    solution: Optional["Solution"] = None,
    processes: Optional[int] = None,
    **kwargs,
) -> pd.DataFrame:
    """Knock out each reaction pair from the combinations of two given lists.

    We say 'pair' here but the order order does not matter.

    Parameters
    ----------
    model : cobra.Model
        The metabolic model to perform deletions in.
    reaction_list1 : list of cobra.Reaction or str, optional
        The first reaction list to be deleted. If not passed,
        all the reactions from the model are used (default None).
    reaction_list2 : list of cobra.Reaction or str, optional
        The second reaction list to be deleted. If not passed,
        all the reactions from the model are used (default None).
    method: {"fba", "moma", "linear moma", "room", "linear room"}, optional
        Method used to predict the growth rate (default "fba").
    solution : cobra.Solution, optional
        A previous solution to use as a reference for (linear) MOMA or ROOM
        (default None).
    processes : int, optional
        The number of parallel processes to run. Can speed up the computations
        if the number of knockouts to perform is large. If not passed,
        will be set to `configuration.processes` (default None).
    **kwargs :
        Keyword arguments are passed on to underlying simulation functions
        such as `add_room`.

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
        **kwargs,
    )


def double_gene_deletion(
    model: Model,
    gene_list1: Optional[List[Union[Gene, str]]] = None,
    gene_list2: Optional[List[Union[Gene, str]]] = None,
    method: str = "fba",
    solution: Optional["Solution"] = None,
    processes: Optional[int] = None,
    **kwargs,
) -> pd.DataFrame:
    """Knock out each gene pair from the combination of two given lists.

    We say 'pair' here but the order order does not matter.

    Parameters
    ----------
    model : cobra.Model
        The metabolic model to perform deletions in.
    gene_list1 : list of cobra.Gene or str, optional
        The first gene list to be deleted. If not passed,
        all the genes from the model are used (default None).
    gene_list2 : list of cobra.Gene or str, optional
        The second gene list to be deleted. If not passed,
        all the genes from the model are used (default None).
    method : {"fba", "moma", "linear moma", "room", "linear room"}, optional
        Method used to predict the growth rate (default None).
    solution : cobra.Solution, optional
        A previous solution to use as a reference for (linear) MOMA or ROOM
        (default None).
    processes : int, optional
        The number of parallel processes to run. Can speed up the computations
        if the number of knockouts to perform is large. If not passed,
        will be set to `configuration.processes` (default None).
    **kwargs :
        Keyword arguments are passed on to underlying simulation functions
        such as `add_room`.

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
        **kwargs,
    )


@pd.api.extensions.register_dataframe_accessor("knockout")
class KnockoutAccessor:
    """
    Access unique combinations of reactions in deletion results.

    This allows acces in the form of `results.knockout[rxn1]` or
    `results.knockout["rxn1_id"]`. Each individual entry will return a
    deletion so `results.knockout[rxn1, rxn2]` will return two deletions
    (for individual knockouts of rxn1 and rxn2 respectively).
    Multi-deletions can be accessed by passing in sets like
    `results.knockout[{rxn1, rxn2}]` which denotes the double deletion of
    both reactions. Thus, the following are allowed index elements:

    - single reactions or genes (depending on whether it is a gene or
      reaction deletion)
    - single reaction IDs or gene IDs
    - lists of single single reaction IDs or gene IDs (will return one row
      for each
      element in the list)
    - sets of reactions or genes (for multi-deletions)
    - sets of reactions IDs or gene IDs
    - list of sets of objects or IDs (to get several multi-deletions)

    Parameters:
    -----------
    pandas_obj : pandas.DataFrame or pandas.Series
        A result from one of the deletion methods.

    """

    def __init__(self, pandas_obj: Union[pd.DataFrame, pd.Series]) -> None:
        """Set up the accessor."""
        self._validate(pandas_obj)
        self._result = pandas_obj

    @staticmethod
    def _validate(obj: pd.DataFrame) -> None:
        """Validate the object given.

        Parameters
        ----------
        obj : pandas.DataFrame
            The object to validate.

        Raises
        ------
        AttributeError
            If the object supplied is not a DataFrame.

        """
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

        Parameters
        ----------
        args : cobra.Reaction, cobra.Gene, str, set, or list
            The deletions to be returned. Accepts:
            - single reactions or genes
            - single reaction IDs or gene IDs
            - lists of single single reaction IDs or gene IDs
            - sets of reactions or genes
            - sets of reactions IDs or gene IDs
            - list of sets of objects or IDs
            See the docs for usage examples.

        Returns
        -------
        pandas.DataFrame
            The deletion result where the chosen entities have been deleted.
            Each row denotes a deletion.

        Raises
        ------
        ValueError
            If any other object is used as index for lookup.

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
                args = [{elem.id for elem in obj} for obj in args]
            except AttributeError:
                args = [set(obj) for obj in args]
        else:
            raise ValueError(
                "Allowed indices are single cobra.Reaction or cobra.Gene, "
                "lists of cobra.Reaction of cobra.Gene, or lists of sets "
                "of cobra.Reaction or cobra.Gene."
            )
        found = [x in args for x in self._result.ids]
        return self._result[found]
