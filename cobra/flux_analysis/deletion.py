# -*- coding: utf-8 -*-

import sys
import multiprocessing
import logging
import optlang
import math
from six import iteritems
from warnings import warn
from itertools import product
from collections import defaultdict
from functools import partial
from builtins import (map, dict)
from future.utils import raise_
from ..manipulation.delete import (find_gene_knockout_reactions)
import cobra.util.solver as sutil


try:
    import scipy
except ImportError:
    moma = None
else:
    from . import moma

LOGGER = logging.getLogger(__name__)


def infinite_defaultdict():
    return defaultdict(infinite_defaultdict)


def defaultdict_to_dict(def_dict):
    if not isinstance(def_dict, defaultdict):
        return def_dict
    for k, v in iteritems(def_dict):
        def_dict[k] = defaultdict_to_dict(v)
    return dict(def_dict)


def biomass_reaction(model):
    return model.objective.expression.as_coefficients_dict()


def restore_biomass(model):
    result = 0
    for k, v in iteritems(model.notes['biomass_reaction']):
        try:
            result += model.reactions.get_by_id(k.name).flux * v
        except KeyError:
            pass
    return result


def _reactions_knockouts_with_restore(model, reactions):
    with model:
        for reaction in reactions:
            reaction.knock_out()
        growth = _get_growth(model)
    return (growth, [r.id for r in reactions])


def _get_growth(model):
    try:
        if 'moma_old_objective' in model.solver.variables:
            model.slim_optimize()
            growth = restore_biomass(model)
        else:
            growth = model.slim_optimize()
        if math.isnan(growth):
            growth = None
    except optlang.exceptions.SolverError:
        growth = None
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
    growth, _ = _reactions_knockouts_with_restore(model, all_reactions)
    return (growth, ids)


def _reaction_deletion_worker(ids):
    global _model
    return _reaction_deletion(_model, ids)


def _gene_deletion_worker(ids):
    global _model
    return _gene_deletion(_model, ids)


def _init_worker(model):
    global _model
    _model = model


def _multi_deletion(cobra_model, entity, element_lists, method="fba",
                    number_of_processes=None, solver=None,
                    zero_cutoff=1e-12):
    """
    Helper function that provides the common interface for sequential
    knockouts

    Parameters
    ----------
    cobra_model : cobra.Model
        The metabolic model to perform deletions in.

    element_lists : list
        List of lists `cobra.Reaction`s or `cobra.Gene`s (or their IDs)
        to be deleted.

    method: {"fba", "moma"}, optional
        Procedure used to predict the growth rate.

    solver: str, optional
        This must be a QP-capable solver for MOMA. If left unspecified,
        a suitable solver will be automatically chosen.

    zero_cutoff: float, optional
        When checking to see if a value is 0, this threshold is used.

    Returns
    -------
    dict of dict
        A sparse representation of all combinations of pairs of reaction
        deletions without replacements.
    """
    with cobra_model as model:
        try:
            (legacy, solver) = sutil.choose_solver(model, solver,
                                                   qp=(method == "moma"))
        except sutil.SolverNotFound:
            if method == "moma":
                warn(
                    "Cannot use MOMA since no QP-capable solver was found. "
                    "Falling back to FBA.")
                (legacy, solver) = sutil.choose_solver(model, solver)
            else:
                (err_type, err_val, err_tb) = sys.exc_info()
                raise_(err_type, err_val, err_tb)  # reraise for Python2&3
        if legacy:
            raise ValueError(
                "Legacy solvers are not supported any longer. "
                "Please use one of the optlang solver interfaces instead.")

        if number_of_processes is None:
            try:
                num_cpu = multiprocessing.cpu_count()
            except NotImplementedError:
                warn("Number of cores could not be detected - assuming 1.")
                num_cpu = 1
        else:
            num_cpu = number_of_processes

        if 'moma' in method:
            model.notes['biomass_reaction'] = biomass_reaction(model)
            moma.add_moma(model, linear='linear' in method)

        args = set([frozenset(comb) for comb in product(*element_lists)])

        def extract_knockout_results(results):
            return {frozenset(ids): 0.0 if growth and abs(
                growth) < zero_cutoff else growth
                    for (growth, ids) in results}

        if num_cpu > 1:
            WORKER_FUNCTIONS = dict(
                gene=_gene_deletion_worker,
                reaction=_reaction_deletion_worker
            )
            chunk_size = len(args) // num_cpu
            pool = multiprocessing.Pool(
                num_cpu, initializer=_init_worker, initargs=(model,)
            )
            results = extract_knockout_results(
                pool.imap_unordered(WORKER_FUNCTIONS[entity], args,
                                    chunksize=chunk_size))
            pool.close()
        else:
            WORKER_FUNCTIONS = dict(
                gene=_gene_deletion,
                reaction=_reaction_deletion
            )
            results = extract_knockout_results(map(
                partial(WORKER_FUNCTIONS[entity], model), args
            ))
        double = infinite_defaultdict()
        for comb in product(*element_lists):
            growth = results[frozenset(comb)]
            current = double
            for v in comb[:-1]:
                current = current[v]
            current[comb[-1]] = growth
        return defaultdict_to_dict(double)


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


def single_reaction_deletion(model, reaction_list=None, **kwargs):
    """Sequentially knocks out each reaction from a given reaction list.

    Parameters
    ----------
    model : a cobra model
        The model from which to delete the reactions. The model will not be
        modified.
    reaction_list : iterable
        List of reaction IDs or cobra.Reaction. If None (default) will use all
        reactions in the model.
    method : str, optional
        The method used to obtain fluxes. Must be one of "fba" or "moma".
    solver : str, optional
        Name of the solver to be used.
    solver_args : optional
        Additional arguments for the solver. Ignored for optlang solver, please
        use `model.solver.configuration` instead.

    Returns
    -------
    tuple of 2 dictionaries
        The first dictionary maps each reaction id to its growth rate after
        the knockout. The second tuple reports the solutions status (for
        instance "optimal" for each knockout).
    """
    return _multi_deletion(
        model,
        'reaction',
        element_lists=_element_lists(model.reactions, reaction_list),
        **kwargs
    )


def single_gene_deletion(model, gene_list=None, **kwargs):
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
    tuple of 2 dictionaries
        The first dictionary maps each gene id to its growth rate after
        the knockout. The second tuple reports the solutions status (for
        instance "optimal" for each knockout).
    """
    return _multi_deletion(
        model,
        'gene',
        element_lists=_element_lists(model.genes, gene_list),
        **kwargs
    )


def double_reaction_deletion(model,
                             reaction_list1=None, reaction_list2=None,
                             **kwargs):
    """
    Sequentially delete pairs of reactions and record the objective value.

    Parameters
    ----------
    model : cobra.Model
        The metabolic model to perform double deletions in.

    reaction_list1, reaction_list2 : list, optional
        Lists of `cobra.Reaction`s (or their IDs) to be deleted. If not
        provided, all reactions in the model sorted by ID will be used.
        `reaction_list2` will use `reaction_list1` by default.

    method: {"fba", "moma"}, optional
        Procedure used to predict the growth rate.

    solver: str, optional
        This must be a QP-capable solver for MOMA. If left unspecified,
        a suitable solver will be automatically chosen.

    zero_cutoff: float, optional
        When checking to see if a value is 0, this threshold is used.

    Returns
    -------
    dict of dict
        A sparse representation of all combinations of pairs of reaction
        deletions without replacements.
    """

    reaction_list1, reaction_list2 = _element_lists(model.reactions,
                                                    reaction_list1,
                                                    reaction_list2)
    return _multi_deletion(model, 'reaction',
                           element_lists=[reaction_list1, reaction_list2],
                           **kwargs)


def double_gene_deletion(model, gene_list1=None, gene_list2=None, **kwargs):
    """Sequentially knocks out pairs of genes in a model

    model : :class:`~cobra.core.Model.Model`
        cobra model in which to perform deletions

    gene_list1 : [:class:`~cobra.core.Gene.Gene`:] (or their id's)
        Genes to be deleted. These will be the keys in the resulting dict.
        If not provided, all genes will be used.

    gene_list1 : [:class:`~cobra.core.Gene.Gene`:] (or their id's)
        Genes to be deleted. These will be the keys in the resulting dict.
        If not provided, gene_list1 will be used.

    method: "fba" or "moma"
        Procedure used to predict the growth rate

    solver: str for solver name
        This must be a QP-capable solver for MOMA. If left unspecified,
        a suitable solver will be automatically chosen.

    zero_cutoff: float
        When checking to see if a value is 0, this threshold is used.

    number_of_processes: int for number of processes to use.
        If unspecified, the number of parallel processes to use will be
        automatically determined. Setting this to 1 explicitly disables used
        of the multiprocessing library.

    .. note:: multiprocessing is not supported with method=moma

    return_frame: bool
        If true, formats the results as a pandas.Dataframe. Otherwise
        returns a dict of the form:
        {"x": row_labels, "y": column_labels", "data": 2D matrix}
    """

    gene_list1, gene_list2 = _element_lists(model.genes, gene_list1,
                                            gene_list2)
    return _multi_deletion(model, 'gene',
                           element_lists=[gene_list1, gene_list2], **kwargs)


def double_deletion(cobra_model, element_list_1=None, element_list_2=None,
                    element_type='gene', **kwargs):
    """Wrapper for double_gene_deletion and double_reaction_deletion

    .. deprecated :: 0.4
        Use double_reaction_deletion and double_gene_deletion
    """
    warn(
        "deprecated - use single_reaction_deletion and single_gene_deletion")
    if element_type == "reaction":
        return double_reaction_deletion(cobra_model, element_list_1,
                                        element_list_2, **kwargs)
    elif element_type == "gene":
        return double_gene_deletion(cobra_model, element_list_1,
                                    element_list_2, **kwargs)
    else:
        raise Exception("unknown element type")
