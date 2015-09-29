import re
from copy import deepcopy
from warnings import warn
from ast import NodeTransformer, And

from six import iteritems, string_types

from ..core.Gene import eval_gpr, parse_gpr, ast2str


def prune_unused_metabolites(cobra_model):
    """Removes metabolites that aren't involved in any reactions in the model

    cobra_model: A Model object.

    """
    inactive_metabolites = []
    active_metabolites = []
    for the_metabolite in cobra_model.metabolites:
        if len(the_metabolite._reaction) == 0:
            the_metabolite.remove_from_model(cobra_model)
            inactive_metabolites.append(the_metabolite)
        else:
            active_metabolites.append(the_metabolite)
    if inactive_metabolites:
        return inactive_metabolites
    else:
        warn('All metabolites used in at least 1 reaction')


def prune_unused_reactions(cobra_model):
    """Removes reactions from cobra_model.

    cobra_model: A Model object.

    reactions_to_prune: None, a string matching a reaction.id, a
    cobra.Reaction, or as list of the ids / Reactions to remove from
    cobra_model.  If None then the function will delete reactions that have no
    active metabolites in the model.

    """
    pruned_reactions = []
    reactions_to_prune = [x for x in cobra_model.reactions
                          if len(x._metabolites) == 0]
    for the_reaction in reactions_to_prune:
        try:
            the_reaction.remove_from_model(cobra_model)
            pruned_reactions.append(the_reaction)
        except:
            warn('%s not in %s' % (the_reaction.id, cobra_model.id))
    if not pruned_reactions:
        warn('All reactions have at least 1 metabolite')
        return


def undelete_model_genes(cobra_model):
    """Undoes the effects of a call to delete_model_genes in place.

    cobra_model:  A cobra.Model which will be modified in place

    """

    if cobra_model._trimmed_genes is not None:
        for x in cobra_model._trimmed_genes:
            x.functional = True

    if cobra_model._trimmed_reactions is not None:
        for the_reaction, (lower_bound, upper_bound) in \
                cobra_model._trimmed_reactions.items():
            the_reaction.lower_bound = lower_bound
            the_reaction.upper_bound = upper_bound

    cobra_model._trimmed_genes = []
    cobra_model._trimmed_reactions = {}
    cobra_model._trimmed = False


def get_compiled_gene_reaction_rules(cobra_model):
    """Generates a dict of compiled gene_reaction_rules

    Any gene_reaction_rule expressions which cannot be compiled or do not
    evaluate after compiling will be excluded. The result can be used in the
    find_gene_knockout_reactions function to speed up evaluation of these
    rules.

    """
    return {r: parse_gpr(r.gene_reaction_rule)[0]
            for r in cobra_model.reactions}


def find_gene_knockout_reactions(cobra_model, gene_list,
                                 compiled_gene_reaction_rules=None):
    """identify reactions which will be disabled when the genes are knocked out

    cobra_model: :class:`~cobra.core.Model.Model`

    gene_list: iterable of :class:`~cobra.core.Gene.Gene`

    compiled_gene_reaction_rules: dict of {reaction_id: compiled_string}
        If provided, this gives pre-compiled gene_reaction_rule strings.
        The compiled rule strings can be evaluated much faster. If a rule
        is not provided, the regular expression evaluation will be used.
        Because not all gene_reaction_rule strings can be evaluated, this
        dict must exclude any rules which can not be used with eval.

    """
    potential_reactions = set()
    for gene in gene_list:
        if isinstance(gene, string_types):
            gene = cobra_model.genes.get_by_id(gene)
        potential_reactions.update(gene._reaction)
    gene_set = {str(i) for i in gene_list}
    if compiled_gene_reaction_rules is None:
        compiled_gene_reaction_rules = {r: parse_gpr(r.gene_reaction_rule)[0]
                                        for r in potential_reactions}

    return [r for r in potential_reactions
            if not eval_gpr(compiled_gene_reaction_rules[r], gene_set)]


def delete_model_genes(cobra_model, gene_list,
                       cumulative_deletions=True, disable_orphans=False):
    """delete_model_genes will set the upper and lower bounds for reactions
    catalysed by the genes in gene_list if deleting the genes means that
    the reaction cannot proceed according to
    cobra_model.reactions[:].gene_reaction_rule

    cumulative_deletions: False or True.  If True then any previous
    deletions will be maintained in the model.

    """
    if disable_orphans:
        raise NotImplementedError("disable_orphans not implemented")
    if not hasattr(cobra_model, '_trimmed'):
        cobra_model._trimmed = False
        cobra_model._trimmed_genes = []
        cobra_model._trimmed_reactions = {}  # Store the old bounds in here.
    # older models have this
    if cobra_model._trimmed_genes is None:
        cobra_model._trimmed_genes = []
    if cobra_model._trimmed_reactions is None:
        cobra_model._trimmed_reactions = {}
    # Allow a single gene to be fed in as a string instead of a list.
    if not hasattr(gene_list, '__iter__') or \
            hasattr(gene_list, 'id'):  # cobra.Gene has __iter__
        gene_list = [gene_list]

    if not hasattr(gene_list[0], 'id'):
        if gene_list[0] in cobra_model.genes:
                tmp_gene_dict = dict([(x.id, x) for x in cobra_model.genes])
        else:
            # assume we're dealing with names if no match to an id
            tmp_gene_dict = dict([(x.name, x) for x in cobra_model.genes])
        gene_list = [tmp_gene_dict[x] for x in gene_list]

    # Make the genes non-functional
    for x in gene_list:
        x.functional = False

    if cumulative_deletions:
        gene_list.extend(cobra_model._trimmed_genes)
    else:
        undelete_model_genes(cobra_model)

    for the_reaction in find_gene_knockout_reactions(cobra_model, gene_list):
        # Running this on an already deleted reaction will overwrite the
        # stored reaction bounds.
        if the_reaction in cobra_model._trimmed_reactions:
            continue
        old_lower_bound = the_reaction.lower_bound
        old_upper_bound = the_reaction.upper_bound
        cobra_model._trimmed_reactions[the_reaction] = (old_lower_bound,
                                                        old_upper_bound)
        the_reaction.lower_bound = 0.
        the_reaction.upper_bound = 0.
        cobra_model._trimmed = True

    cobra_model._trimmed_genes = list(set(cobra_model._trimmed_genes +
                                          gene_list))


class _GeneRemover(NodeTransformer):
    def __init__(self, target_genes):
        NodeTransformer.__init__(self)
        self.target_genes = {str(i) for i in target_genes}

    def visit_Name(self, node):
        return None if node.id in self.target_genes else node

    def visit_BoolOp(self, node):
        original_n = len(node.values)
        self.generic_visit(node)
        if len(node.values) == 0:
            return None
        # AND with any entities removed
        if len(node.values) < original_n and isinstance(node.op, And):
            return None
        # if one entity in an OR was removed, just that entity passed up
        if len(node.values) == 1:
            return node.values[0]
        return node


def remove_genes(cobra_model, gene_list, remove_reactions=True):
    """remove genes entirely from the model

    This will also simplify all gene_reaction_rules with this
    gene inactivated."""
    gene_set = {cobra_model.genes.get_by_id(str(i)) for i in gene_list}
    gene_id_set = {i.id for i in gene_set}
    remover = _GeneRemover(gene_id_set)
    ast_rules = get_compiled_gene_reaction_rules(cobra_model)
    target_reactions = []
    for reaction, rule in iteritems(ast_rules):
        if reaction.gene_reaction_rule is None or \
                len(reaction.gene_reaction_rule) == 0:
            continue
        # reactions to remove
        if remove_reactions and not eval_gpr(rule, gene_id_set):
            target_reactions.append(reaction)
        else:
            # if the reaction is not removed, remove the gene
            # from its gpr
            remover.visit(rule)
            new_rule = ast2str(rule)
            if new_rule != reaction.gene_reaction_rule:
                reaction.gene_reaction_rule = new_rule
    for gene in gene_set:
        cobra_model.genes.remove(gene)
    cobra_model.remove_reactions(target_reactions)
