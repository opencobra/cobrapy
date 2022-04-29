"""Provide functions for pruning reactions, metabolites and genes."""
import logging
from ast import And, BoolOp, Name, NodeTransformer
from functools import partial
from typing import TYPE_CHECKING, Iterable, List, Optional, Set, Tuple, Union

from cobra.util import get_context


if TYPE_CHECKING:
    from cobra import Gene, Metabolite, Model, Reaction

logger = logging.getLogger(__name__)


def prune_unused_metabolites(model: "Model") -> Tuple["Model", List["Metabolite"]]:
    """Remove metabolites not involved in any reactions.

    Parameters
    ----------
    model: cobra.Model
        The model to remove unused metabolites from.

    Returns
    -------
    cobra.Model
        Input model with unused metabolites removed.
    list of cobra.Metabolite
        List of metabolites that were removed.

    """
    output_model = model.copy()
    inactive_metabolites = [
        m for m in output_model.metabolites if len(m.reactions) == 0
    ]
    output_model.remove_metabolites(inactive_metabolites)
    return output_model, inactive_metabolites


def prune_unused_reactions(model: "Model") -> Tuple["Model", List["Reaction"]]:
    """Remove reactions with no assigned metabolites, returns pruned model.

    Parameters
    ----------
    model: cobra.Model
        The model to remove unused reactions from.

    Returns
    -------
    cobra.Model
        Input model with unused reactions removed.
    list of cobra.Reaction
        List of reactions that were removed.

    """
    output_model = model.copy()
    reactions_to_prune = [r for r in output_model.reactions if len(r.metabolites) == 0]
    output_model.remove_reactions(reactions_to_prune)
    return output_model, reactions_to_prune


def knock_out_model_genes(
    model: "Model",
    gene_list: Iterable[Union["Gene", int, str]],
) -> List["Reaction"]:
    """Disable the genes in `gene_list`.

    It sets the bounds to "zero" for reactions catalysed by the genes in
    `gene_list` if deleting the genes would stop the reactions from proceeding.

    The changes are reverted upon exit if executed within the model as context.

    Parameters
    ----------
    model: cobra.Model
        The model whose reaction bounds are to be set.
    gene_list: list of cobra.Gene
        The list of genes to knock-out.

    Returns
    -------
    list[cobra.Reaction]
        A list of cobra.Reactions that had the bounds turned to zero.
    """
    rxn_set = set()
    for gene in model.genes.get_by_any(gene_list):
        gene.knock_out()
        rxn_set.update(gene.reactions)
    return [rxn for rxn in rxn_set if not rxn.functional]


def delete_model_genes(
    model: "Model",
    gene_list: Union[List["Gene"], Set["Gene"], List[str], Set[str]],
    cumulative_deletions: bool = True,
    disable_orphans: bool = False,
) -> None:
    """Temporarily remove the effect of genes in `gene_list`.

    It sets the bounds to "zero" for reactions catalysed by the genes in
    `gene_list` if deleting the genes stops the reactions from proceeding.

    Parameters
    ----------
    model: cobra.Model
        The model whose reaction bounds are to be set.
    gene_list: list of cobra.Gene
        The list of genes to knock-out.
    cumulative_deletions: bool, optional
        If True, then any previous deletions will be maintained in the
        model (default True). Unused, ignored.
    disable_orphans: bool, optional
        If True, then orphan reactions will be disabled. Currently, this
        is not implemented (default False). Unused, ignored.

    .. deprecated :: 0.25
        Use cobra.manipulation.knock_out_model_genes to simulate knockouts
        and cobra.manipulation.remove_genes to remove genes from
        the model.

    See Also
    --------
    knock_out_model_genes()
    remove_model_genes()

    """
    logger.warning(
        "delete_model_genes and undelete_model_genes were deprecated "
        "and should not be used. undelete_model_genes was removed."
        "Use cobra.manipulation.remove_genes instead to remove genes "
        "from the model."
    )
    logger.warning(
        f"Use cobra.manipulation.knock_out_model_genes({model}, {gene_list})"
        f" to cause knockouts."
    )
    logger.warning(
        f"If you want to temporarily knock out genes (like undelete_model_genes) did, "
        f"run it in model context, for example"
        f"with({model}):"
        f"    knocked_out_rxns = knock_out_model_genes({model}, {gene_list})"
        f" since it will be reset when leaving the context."
    )


class _GeneRemover(NodeTransformer):
    """
    Class to represent a gene set remover.

    Parameters
    ----------
    target_genes: list or set of cobra.Gene
        A set of genes to be removed.

    """

    def __init__(self, target_genes: Set["Gene"], **kwargs) -> None:
        """Initialize a new object.

        Other Parameters
        ----------------
        kwargs:
            Further keyword arguments are passed on to the parent class.

        """
        super().__init__(**kwargs)
        self.target_genes = {str(i) for i in target_genes}

    def visit_Name(self, node: "Name") -> Optional["Name"]:
        """Remove a gene.

        Parameters
        ----------
        node: ast.Name
            The gene to remove.

        Returns
        -------
        cobra.Gene or None
            None if gene object is in `target_genes`.

        """
        return None if node.id in self.target_genes else node

    def visit_BoolOp(self, node: "BoolOp") -> Optional[Union["BoolOp", "Name"]]:
        """Rules for boolean operations.

        Parameters
        ----------
        node: ast.Name
            The gene to apply rules to.

        Returns
        -------
        ast.Name or None
            None if size of Or node values is zero after applying rule,
            or size of And node values is lower after applying rule.

        """
        original_n = len(node.values)
        self.generic_visit(node)
        if len(node.values) == 0:
            return None
        # AND with any entities removed
        if len(node.values) < original_n and isinstance(node.op, And):
            return None
        # if one entity in an OR was left, just that entity passed up
        if len(node.values) == 1:
            return node.values[0]
        return node


def remove_genes(
    model: "Model",
    gene_list: Union[List["Gene"], Set["Gene"], List[str], Union[str]],
    remove_reactions: bool = True,
) -> None:
    """Remove genes entirely from the model.

    This will also simplify all gene-reaction rules with the genes
    inactivated.

    Parameters
    ----------
    model: cobra.Model
        The model to remove genes from.
    gene_list: list of cobra.Gene or gene ids
        The list of gene objects to remove.
    remove_reactions: bool, optional
        Whether to remove reactions associated with genes in `gene_list`
        (default True).

    """
    gene_set = {model.genes.get_by_id(str(i)) for i in gene_list}
    gene_id_set = {i.id for i in gene_set}
    remover = _GeneRemover(gene_id_set)
    target_reactions = []
    rxns_to_revisit = set()
    context = get_context(model)
    for rxn in model.reactions:
        if rxn.gene_reaction_rule is None or len(rxn.gene_reaction_rule) == 0:
            continue
        # reactions to remove
        if remove_reactions and not rxn.gpr.eval(gene_id_set):
            target_reactions.append(rxn)
        else:
            # if the reaction is not removed, remove the gene
            # from its gpr
            old_gpr = rxn.gpr.copy()
            remover.visit(rxn.gpr)
            # If the remover completely removed the AST tree from the GPR, it will not
            # have body at all, which is why this isn't if body is None.
            if not hasattr(rxn.gpr, "body"):
                rxn.gpr.body = None
                rxn._genes = set()
            else:
                rxns_to_revisit.add(rxn)
            if context:
                context(partial(setattr, rxn, "gpr", old_gpr))
                context(partial(rxn.update_genes_from_gpr))
    for gene in gene_set:
        model.genes.remove(gene)
        gene._model = None
        if context:
            context(partial(model.genes.add, gene))
            context(partial(setattr, gene, "_model", model))
        # remove reference to the gene in all groups
        associated_groups = model.get_associated_groups(gene)
        for group in associated_groups:
            group.remove_members(gene)
    model.remove_reactions(target_reactions)
    for rxn in rxns_to_revisit:
        rxn.update_genes_from_gpr()
