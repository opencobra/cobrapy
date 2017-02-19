# -*- coding: utf-8 -*-

import sys
import multiprocessing
import logging
from warnings import warn
from itertools import (chain, product)
from operator import attrgetter
from collections import defaultdict
from builtins import (map, dict, zip)
from future.utils import raise_

from six import iteritems
import numpy
from tqdm import tqdm

from ..solvers import get_solver_name, solver_dict
from ..manipulation.delete import (find_gene_knockout_reactions,
                                   get_compiled_gene_reaction_rules)
from .deletion_worker import CobraDeletionPool, CobraDeletionMockPool
import cobra.util.solver as sutil

try:
    import scipy
except ImportError:
    moma = None
else:
    from . import moma

try:
    from pandas import DataFrame
except:
    DataFrame = None

LOGGER = logging.getLogger(__name__)


def double_reaction_deletion(model,
                             reaction_list1=None, reaction_list2=None,
                             method="fba", num_cpu=None,
                             solver=None, zero_cutoff=1e-12,
                             **kwargs):
    """
    Sequentially delete pairs of reactions and record the objective value.

    Parameters
    ----------
    cobra_model : cobra.Model
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
    try:
        (legacy, solver) = sutil.choose_solver(model, solver, qp=(method == "moma"))
    except sutil.SolverNotFound:
        if method == "moma":
            warn("Cannot use MOMA since no QP-capable solver was found. Falling"
                 " back to FBA.")
            (legacy, solver) = sutil.choose_solver(model, solver)
        else:
            (err_type, err_val, err_tb) = sys.exc_info()
            raise_(err_type, err_val, err_tb)  # reraise for Python2&3

    if legacy:
        raise ValueError(
            "Legacy solvers are not supported any longer. Please use one of"
            " the optlang solver interfaces instead.")

    if num_cpu is None:
        num_cpu = kwargs.get("number_of_processes")
    if num_cpu is None:
        try:
            num_cpu = multiprocessing.cpu_count()
        except NotImplementedError:
            warn("Number of cores could not be detected - assuming 1.")
            num_cpu = 1

    print("making reaction lists")

    if reaction_list1 is None:
        reaction_list1 = sorted(model.reactions, key=attrgetter("id"))
    else:
        reaction_list1 = model.reactions.get_by_any(reaction_list1)

    if reaction_list2 is None:
        reaction_list2 = reaction_list1
    else:
        reaction_list2 = model.reactions.get_by_any(reaction_list2)

    print("making first solution")
    # Reactions with 0 flux in the current conditions have no impact on
    # deletion. Determine those reactions now to speed up computation later.
    solution = model.optimize()
    if solution.status == "optimal":
        kwargs["growth"] = solution.objective_value
        kwargs["zero_flux_reactions"] = \
            {rxn for (rxn, flux) in solution.fluxes.items()
             if abs(flux) < zero_cutoff}
    else:
        warn("Non-optimal ({0:s}) base solution.".format(solution.status))

    print("making args")
    args = list()
    covered = set()
#    for rxn_a in reaction_list1:
#        id_a = rxn_a.id
#        for rxn_b in reaction_list2:
#            if rxn_a is rxn_b:
#                continue
#            id_b = rxn_b.id
#            if ((id_a, id_b) in covered) or ((id_b, id_a) in covered):
#                continue
#            args.append((model, id_a, id_b))
#            covered.add((id_a, id_b))

    for rxn_a in reaction_list1:
        id_a = rxn_a.id
        rest = list()
        for rxn_b in reaction_list2:
            if rxn_a is rxn_b:
                continue
            id_b = rxn_b.id
            if ((id_a, id_b) in covered) or ((id_b, id_a) in covered):
                continue
            rest.append(id_b)
            covered.add((id_a, id_b))
        args.append((model, id_a, rest))

    print("{0:n} double knock-outs".format(len(covered)))

    print("making iterators")
    if num_cpu > 1:
        # multiprocessing
#        chunk_size = len(covered) // (num_cpu * 2)
        chunk_size = 1
        pool = multiprocessing.Pool(num_cpu)
        results = pool.imap_unordered(_deletion_worker2, args,
            chunksize=chunk_size)
    else:
        results = map(_deletion_worker2, args)

    print("saving results")
    double = defaultdict(dict)

#    for (growth, id_a, id_b) in results:
#        if growth is None:
#            double[id_a][id_b] = None
#        continue
#        double[id_a][id_b] = 0.0 if abs(growth) < zero_cutoff else growth

    with tqdm(total=len(covered)) as pbar:
        for (growth, prim_id, rest) in results:
            for (rxn_id, value) in zip(rest, growth):
                pbar.update()
                if value is None:
                    double[prim_id][rxn_id] = None
                continue
                double[prim_id][rxn_id] = 0.0 if abs(value) < zero_cutoff else value

    return dict(double)

def _deletion_worker(args):
    model = args[0]
    rxn_ids = args[1:]
    with model as mod:
        for rxn_id in rxn_ids:
            mod.reactions.get_by_id(rxn_id).knock_out()
        try:
            solution = mod.optimize()
            growth = solution.objective_value
        except:
#            LOGGER.error(mod.solver.status)
            growth = None
    return (growth,) + tuple(rxn_ids)

def _deletion_worker2(args):
    (model, primary_id, rest) = args
    primary = model.reactions.get_by_id(primary_id)
    growth = list()
    with model as mod:
        primary.knock_out()
        for rxn_id in rest:
            with mod as double:
                double.reactions.get_by_id(rxn_id).knock_out()
                try:
                    solution = double.optimize()
                    growth.append(solution.objective_value)
                except:
#                    LOGGER.error(double.solver.status)
                    growth.append(None)
    return (growth, primary_id, rest)


# Utility functions
def generate_matrix_indexes(ids1, ids2):
    """map an identifier to an entry in the square result matrix"""
    return {id: index for index, id in enumerate(set(chain(ids1, ids2)))}


def yield_upper_tria_indexes(ids1, ids2, id_to_index):
    """gives the necessary indexes in the upper triangle

    ids1 and ids2 are lists of the identifiers i.e. gene id's or reaction
    indexes to be knocked out. id_to_index maps each identifier to its index
    in the result matrix.

    Note that this does not return indexes for the diagonal. Those have
    to be computed separately."""
    # sets to check for inclusion in o(1)
    id_set1 = set(ids1)
    id_set2 = set(ids2)
    for id1, id2 in product(ids1, ids2):
        # indexes in the result matrix
        index1 = id_to_index[id1]
        index2 = id_to_index[id2]
        # upper triangle
        if index2 > index1:
            yield ((index1, index2), (id1, id2))
        # lower triangle but would be skipped, so return in upper triangle
        elif id2 not in id_set1 or id1 not in id_set2:
            yield((index2, index1), (id2, id1))  # note that order flipped


def _format_upper_triangular_matrix(row_indexes, column_indexes, matrix):
    """reformat the square upper-triangular result matrix

    For example, results may look like this
    [[ A  B  C  D]
     [ -  -  -  -]
     [ -  -  E  F]
     [ -  -  -  G]]
    In this case, the second row was skipped. This means we have
    row_indexes [0, 2, 3] and column_indexes [0, 1, 2, 3]

    First, it will reflect the upper triangle into the lower triangle
    [[ A  B  C  D]
     [ B  -  -  -]
     [ C  -  E  F]
     [ D  -  F  G]]

    Finally, it will remove the missing rows and return
    [[ A  B  C  D]
     [ C  -  E  F]
     [ D  -  F  G]]
    """
    results = matrix.copy()
    # Thse select the indexes for the upper triangle. However, switching
    # the order selects the lower triangle.
    triu1, triu2 = numpy.triu_indices(matrix.shape[0])
    # This makes reflection pretty easy
    results[triu2, triu1] = results[triu1, triu2]
    # Remove the missing rows and  return.
    return results[row_indexes, :][:, column_indexes]


def format_results_frame(row_ids, column_ids, matrix, return_frame=False):
    """format results as a pandas.DataFrame if desired/possible

    Otherwise returns a dict of
    {"x": row_ids, "y": column_ids", "data": result_matrx}"""
    if return_frame and DataFrame:
        return DataFrame(data=matrix, index=row_ids, columns=column_ids)
    elif return_frame and not DataFrame:
        warn("could not import pandas.DataFrame")
    return {"x": row_ids, "y": column_ids, "data": matrix}


def double_deletion(cobra_model, element_list_1=None, element_list_2=None,
                    element_type='gene', **kwargs):
    """Wrapper for double_gene_deletion and double_reaction_deletion

    .. deprecated :: 0.4
        Use double_reaction_deletion and double_gene_deletion
    """
    warn("deprecated - use single_reaction_deletion and single_gene_deletion")
    if element_type == "reaction":
        return double_reaction_deletion(cobra_model, element_list_1,
                                        element_list_2, **kwargs)
    elif element_type == "gene":
        return double_gene_deletion(cobra_model, element_list_1,
                                    element_list_2, **kwargs)
    else:
        raise Exception("unknown element type")


#def double_reaction_deletion(cobra_model,
#                             reaction_list1=None, reaction_list2=None,
#                             method="fba", return_frame=False,
#                             solver=None, zero_cutoff=1e-12,
#                             **kwargs):
#    """sequentially knocks out pairs of reactions in a model
#
#    cobra_model : :class:`~cobra.core.Model.Model`
#        cobra model in which to perform deletions
#
#    reaction_list1 : [:class:`~cobra.core.Reaction.Reaction`:] (or their id's)
#        Reactions to be deleted. These will be the rows in the result.
#        If not provided, all reactions will be used.
#
#    reaction_list2 : [:class:`~cobra.core.Reaction`:] (or their id's)
#        Reactions to be deleted. These will be the rows in the result.
#        If not provided, reaction_list1 will be used.
#
#    method: "fba" or "moma"
#        Procedure used to predict the growth rate
#
#    solver: str for solver name
#        This must be a QP-capable solver for MOMA. If left unspecified,
#        a suitable solver will be automatically chosen.
#
#    zero_cutoff: float
#        When checking to see if a value is 0, this threshold is used.
#
#    return_frame: bool
#        If true, formats the results as a pandas.Dataframe. Otherwise
#        returns a dict of the form:
#        {"x": row_labels, "y": column_labels", "data": 2D matrix}
#    """
#    # handle arguments which need to be passed on
#    if solver is None:
#        solver = get_solver_name(qp=(method == "moma"))
#        kwargs["solver"] = solver
#    kwargs["zero_cutoff"] = zero_cutoff
#
#    # generate other arguments
#
#    # identifiers for reactions are their indexes
#    if reaction_list1 is None:
#        reaction_indexes1 = range(len(cobra_model.reactions))
#    else:
#        reaction_indexes1 = [cobra_model.reactions.index(r)
#                             for r in reaction_list1]
#    if reaction_list2 is None:
#        reaction_indexes2 = reaction_indexes1
#    else:
#        reaction_indexes2 = [cobra_model.reactions.index(r)
#                             for r in reaction_list2]
#    reaction_to_result = generate_matrix_indexes(reaction_indexes1,
#                                                 reaction_indexes2)
#
#    # Determine 0 flux reactions. If an optimal solution passes no flux
#    # through the deleted reactions, then we know removing them will
#    # not change the solution.
#    wt_solution = solver_dict[solver].solve(cobra_model)
#    if wt_solution.status == "optimal":
#        kwargs["wt_growth_rate"] = wt_solution.f
#        kwargs["no_flux_reaction_indexes"] = \
#            {i for i, v in enumerate(wt_solution.x) if abs(v) < zero_cutoff}
#    else:
#        warn("wild-type solution status is '%s'" % wt_solution.status)
#
#    # call the computing functions
#    if method == "fba":
#        results = _double_reaction_deletion_fba(
#            cobra_model, reaction_indexes1, reaction_indexes2,
#            reaction_to_result, **kwargs)
#    elif method == "moma":
#        results = _double_reaction_deletion_moma(
#            cobra_model, reaction_indexes1, reaction_indexes2,
#            reaction_to_result, **kwargs)
#    else:
#        raise ValueError("Unknown deletion method '%s'" % method)
#
#    # convert upper triangular matrix to full matrix
#    full_result = _format_upper_triangular_matrix(
#        [reaction_to_result[i] for i in reaction_indexes1],  # row indexes
#        [reaction_to_result[i] for i in reaction_indexes2],  # col indexes
#        results)
#
#    # format appropriately with labels
#    row_ids = [cobra_model.reactions[i].id for i in reaction_indexes1]
#    column_ids = [cobra_model.reactions[i].id for i in reaction_indexes2]
#    return format_results_frame(row_ids, column_ids,
#                                full_result, return_frame)


def double_gene_deletion(cobra_model,
                         gene_list1=None, gene_list2=None,
                         method="fba", return_frame=False,
                         solver=None, zero_cutoff=1e-12,
                         **kwargs):
    """sequentially knocks out pairs of genes in a model

    cobra_model : :class:`~cobra.core.Model.Model`
        cobra model in which to perform deletions

    gene_list1 : [:class:`~cobra.core.Gene.Gene`:] (or their id's)
        Genes to be deleted. These will be the rows in the result.
        If not provided, all reactions will be used.

    gene_list1 : [:class:`~cobra.core.Gene.Gene`:] (or their id's)
        Genes to be deleted. These will be the rows in the result.
        If not provided, reaction_list1 will be used.

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
    # handle arguments which need to be passed on
    if solver is None:
        solver = get_solver_name(qp=(method == "moma"))
        kwargs["solver"] = solver
    kwargs["zero_cutoff"] = zero_cutoff

    # generate other arguments

    # identifiers for genes
    if gene_list1 is None:
        gene_ids1 = cobra_model.genes.list_attr("id")
    else:
        gene_ids1 = [str(i) for i in gene_list1]
    if gene_list2 is None:
        gene_ids2 = gene_ids1
    else:
        gene_ids2 = [str(i) for i in gene_list2]

    # The gene_id_to_result dict will map each gene id to the index
    # in the result matrix.
    gene_id_to_result = generate_matrix_indexes(gene_ids1, gene_ids2)

    # Determine 0 flux reactions. If an optimal solution passes no flux
    # through the deleted reactions, then we know removing them will
    # not change the solution.
    wt_solution = solver_dict[solver].solve(cobra_model)
    if wt_solution.status == "optimal":
        kwargs["wt_growth_rate"] = wt_solution.f
        kwargs["no_flux_reaction_indexes"] = \
            {i for i, v in enumerate(wt_solution.x) if abs(v) < zero_cutoff}
    else:
        warn("wild-type solution status is '%s'" % wt_solution.status)

    if method == "fba":
        result = _double_gene_deletion_fba(cobra_model, gene_ids1, gene_ids2,
                                           gene_id_to_result, **kwargs)
    elif method == "moma":
        result = _double_gene_deletion_moma(cobra_model, gene_ids1, gene_ids2,
                                            gene_id_to_result, **kwargs)
    else:
        raise ValueError("Unknown deletion method '%s'" % method)

    # convert upper triangular matrix to full matrix
    full_result = _format_upper_triangular_matrix(
        [gene_id_to_result[id] for id in gene_ids1],  # row indexes
        [gene_id_to_result[id] for id in gene_ids2],  # col indexes,
        result)

    # format as a Dataframe if required
    return format_results_frame(gene_ids1, gene_ids2,
                                full_result, return_frame)


def _double_reaction_deletion_fba(cobra_model, reaction_indexes1,
                                  reaction_indexes2, reaction_to_result,
                                  solver, number_of_processes=None,
                                  zero_cutoff=1e-15, wt_growth_rate=None,
                                  no_flux_reaction_indexes=set(), **kwargs):
    """compute double reaction deletions using fba

    cobra_model: model

    reaction_indexes1, reaction_indexes2: reaction indexes (used as unique
        identifiers)

    reaction_to_result: maps each reaction identifier to the entry in
        the result matrix

    no_flux_reaction_indexes: set of indexes for reactions in the model
        which carry no flux in an optimal solution. For deletions only in
        this set, the result will beset to wt_growth_rate.

    returns an upper triangular square matrix
    """
    if solver is None:
        solver = get_solver_name()

    # generate the square result matrix
    n_results = len(reaction_to_result)
    results = numpy.empty((n_results, n_results))
    results.fill(numpy.nan)

    PoolClass = CobraDeletionMockPool if number_of_processes == 1 \
        else CobraDeletionPool  # explicitly disable multiprocessing

    with PoolClass(cobra_model, n_processes=number_of_processes,
                   solver=solver, **kwargs) as pool:

        # precompute all single deletions in the pool and store them along
        # the diagonal
        for reaction_index, result_index in iteritems(reaction_to_result):
            pool.submit((reaction_index, ), label=result_index)
        for result_index, value in pool.receive_all():
            # if singly lethal, set everything in row and column to 0
            value = value if abs(value) > zero_cutoff else 0.
            if value == 0.:
                results[result_index, :] = 0.
                results[:, result_index] = 0.
            else:  # only the diagonal needs to be set
                results[result_index, result_index] = value

        # Run double knockouts in the upper triangle
        index_selector = yield_upper_tria_indexes(
            reaction_indexes1, reaction_indexes2, reaction_to_result)
        for result_index, (r1_index, r2_index) in index_selector:
            # skip if the result was already computed to be lethal
            if results[result_index] == 0:
                continue
            # reactions removed carry no flux
            if r1_index in no_flux_reaction_indexes and \
                    r2_index in no_flux_reaction_indexes:
                results[result_index] = wt_growth_rate
                continue
            pool.submit((r1_index, r2_index), label=result_index)
        # get results
        for result in pool.receive_all():
            results[result[0]] = result[1]

    return results


def _double_gene_deletion_fba(cobra_model, gene_ids1, gene_ids2,
                              gene_id_to_result, solver,
                              number_of_processes=None, zero_cutoff=1e-12,
                              wt_growth_rate=None,
                              no_flux_reaction_indexes=set(), **kwargs):
    """compute double gene deletions using fba

    cobra_model: model

    gene_ids1, gene_ids2: lists of id's to be knocked out

    gene_id_to_result: maps each gene identifier to the entry in
        the result matrix

    no_flux_reaction_indexes: set of indexes for reactions in the model
        which carry no flux in an optimal solution. For deletions only in
        this set, the result will beset to wt_growth_rate.

    returns an upper triangular square matrix
    """
    # Because each gene reaction rule will be evaluated multiple times
    # the reaction has multiple associated genes being deleted, compiling
    # the gene reaction rules ahead of time increases efficiency greatly.
    compiled_rules = get_compiled_gene_reaction_rules(cobra_model)

    n_results = len(gene_id_to_result)
    results = numpy.empty((n_results, n_results))
    results.fill(numpy.nan)

    if number_of_processes == 1:  # explicitly disable multiprocessing
        PoolClass = CobraDeletionMockPool
    else:
        PoolClass = CobraDeletionPool
    with PoolClass(cobra_model, n_processes=number_of_processes,
                   solver=solver, **kwargs) as pool:
        # precompute all single deletions in the pool and store them along
        # the diagonal
        for gene_id, gene_result_index in iteritems(gene_id_to_result):
            ko_reactions = find_gene_knockout_reactions(
                cobra_model, (cobra_model.genes.get_by_id(gene_id),))
            ko_indexes = [cobra_model.reactions.index(i) for i in ko_reactions]
            pool.submit(ko_indexes, label=gene_result_index)
        for result_index, value in pool.receive_all():
            # if singly lethal, set everything in row and column to 0
            value = value if abs(value) > zero_cutoff else 0.
            if value == 0.:
                results[result_index, :] = 0.
                results[:, result_index] = 0.
            else:  # only the diagonal needs to be set
                results[result_index, result_index] = value

        # Run double knockouts in the upper triangle
        index_selector = yield_upper_tria_indexes(gene_ids1, gene_ids2,
                                                  gene_id_to_result)
        for result_index, (gene1, gene2) in index_selector:

            # if singly lethal the results have already been set
            if results[result_index] == 0:
                continue
            ko_reactions = find_gene_knockout_reactions(
                cobra_model, (gene1, gene2), compiled_rules)
            ko_indexes = [cobra_model.reactions.index(i)
                          for i in ko_reactions]
            # if all removed gene indexes carry no flux
            if len(set(ko_indexes) - no_flux_reaction_indexes) == 0:
                results[result_index] = wt_growth_rate
                continue
            pool.submit(ko_indexes, label=result_index)

        for result in pool.receive_all():
            value = result[1]
            if value < zero_cutoff:
                value = 0
            results[result[0]] = value

    return results


def _double_reaction_deletion_moma(cobra_model, reaction_indexes1,
                                   reaction_indexes2, reaction_to_result,
                                   solver, number_of_processes=1,
                                   zero_cutoff=1e-15, wt_growth_rate=None,
                                   no_flux_reaction_indexes=set(), **kwargs):
    """compute double reaction deletions using moma

    cobra_model: model

    reaction_indexes1, reaction_indexes2: reaction indexes (used as unique
        identifiers)

    reaction_to_result: maps each reaction identifier to the entry in
        the result matrix

    no_flux_reaction_indexes: set of indexes for reactions in the model
        which carry no flux in an optimal solution. For deletions only in
        this set, the result will beset to wt_growth_rate.

    number_of_processes: must be 1. Parallel MOMA not yet implmemented

    returns an upper triangular square matrix
    """
    if number_of_processes > 1:
        raise NotImplementedError("parallel MOMA not implemented")
    if moma is None:
        raise RuntimeError("scipy required for MOMA")

    # generate the square result matrix
    n_results = len(reaction_to_result)
    results = numpy.empty((n_results, n_results))
    results.fill(numpy.nan)

    # function to compute reaction knockouts with moma
    moma_model, moma_obj = moma.create_euclidian_moma_model(cobra_model)

    def run(indexes):
        # If all the reactions carry no flux, deletion will have no effect.
        if no_flux_reaction_indexes.issuperset(indexes):
            return wt_growth_rate
        return moma.moma_knockout(moma_model, moma_obj, indexes,
                                  solver=solver, **kwargs).f

    # precompute all single deletions and store them along the diagonal
    for reaction_index, result_index in iteritems(reaction_to_result):
        value = run((reaction_index,))
        value = value if abs(value) > zero_cutoff else 0.
        results[result_index, result_index] = value
        # if singly lethal, the entire row and column are set to 0
        if value == 0.:
            results[result_index, :] = 0.
            results[:, result_index] = 0.

    # Run double knockouts in the upper triangle
    index_selector = yield_upper_tria_indexes(
        reaction_indexes1, reaction_indexes2, reaction_to_result)
    for result_index, (r1_index, r2_index) in index_selector:
        # skip if the result was already computed to be lethal
        if results[result_index] == 0:
            continue
        else:
            results[result_index] = run((r1_index, r2_index))

    return results


def _double_gene_deletion_moma(cobra_model, gene_ids1, gene_ids2,
                               gene_id_to_result, solver,
                               number_of_processes=1,
                               zero_cutoff=1e-12, wt_growth_rate=None,
                               no_flux_reaction_indexes=set(), **kwargs):
    """compute double gene deletions using moma

    cobra_model: model

    gene_ids1, gene_ids2: lists of id's to be knocked out

    gene_id_to_result: maps each gene identifier to the entry in
        the result matrix

    number_of_processes: must be 1. Parallel MOMA not yet implemented

    no_flux_reaction_indexes: set of indexes for reactions in the model
        which carry no flux in an optimal solution. For deletions only in
        this set, the result will beset to wt_growth_rate.

    returns an upper triangular square matrix
    """
    if number_of_processes > 1:
        raise NotImplementedError("parallel MOMA not implemented")
    if moma is None:
        raise RuntimeError("scipy required for MOMA")

    # Because each gene reaction rule will be evaluated multiple times
    # the reaction has multiple associated genes being deleted, compiling
    # the gene reaction rules ahead of time increases efficiency greatly.
    compiled_rules = get_compiled_gene_reaction_rules(cobra_model)

    # function to compute reaction knockouts with moma
    moma_model, moma_obj = moma.create_euclidian_moma_model(cobra_model)

    def run(gene_ids):
        ko_reactions = find_gene_knockout_reactions(cobra_model, gene_ids)
        ko_indexes = map(cobra_model.reactions.index, ko_reactions)
        # If all the reactions carry no flux, deletion will have no effect.
        if no_flux_reaction_indexes.issuperset(gene_ids):
            return wt_growth_rate
        return moma.moma_knockout(moma_model, moma_obj, ko_indexes,
                                  solver=solver, **kwargs).f

    n_results = len(gene_id_to_result)
    results = numpy.empty((n_results, n_results))
    results.fill(numpy.nan)

    # precompute all single deletions and store them along the diagonal
    for gene_id, result_index in iteritems(gene_id_to_result):
        value = run((gene_id,))
        value = value if abs(value) > zero_cutoff else 0.
        results[result_index, result_index] = value
        # If singly lethal, the entire row and column are set to 0.
        if value == 0.:
            results[result_index, :] = 0.
            results[:, result_index] = 0.

    # Run double knockouts in the upper triangle
    index_selector = yield_upper_tria_indexes(gene_ids1, gene_ids2,
                                              gene_id_to_result)
    for result_index, (gene1, gene2) in index_selector:
        # if singly lethal the results have already been set
        if results[result_index] == 0:
            continue
        results[result_index] = run((gene1, gene2))

    return results
