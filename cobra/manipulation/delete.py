import re
from copy import deepcopy
from warnings import warn


# compile regular expressions now instead of in every function call
spontaneous_re = re.compile('(^|(?<=( |\()))s0001(?=( |\)|$))')
and_re = re.compile(r'\band\b')
or_re = re.compile(r'\bor\b')


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
    rules = {}
    # some gene names can not be put through eval
    bad_genes = {g for g in cobra_model.genes if g.id[0].isdigit() or
                 not g.id.isalnum()}
    for reaction in cobra_model.reactions:
        try:
            if len(bad_genes.intersection(reaction.genes)) > 0:
                continue
            rules[reaction.id] = compile(reaction.gene_reaction_rule,
                                         '<string>', 'eval')
        except SyntaxError:
            # This is necessary because some gene_reaction_rules do not
            # compile.
            None
    return rules


def find_gene_knockout_reactions(cobra_model, gene_list,
                                 compiled_gene_reaction_rules={}):
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
    for x in gene_list:
        potential_reactions.update(x._reaction)

    knocked_out_reactions = []
    for the_reaction in potential_reactions:
        # Attempt to use the compiled gene reaction rule if provided
        if the_reaction.id in compiled_gene_reaction_rules:
            gene_state = {i.id: False if i in gene_list else True
                          for i in the_reaction._genes}
            result = eval(compiled_gene_reaction_rules[the_reaction.id],
                          {}, gene_state)
            if result is False:
                knocked_out_reactions.append(the_reaction)
                continue
            elif result is True:
                continue

        # operates on a copy
        gene_reaction_rule = and_re.sub("*", the_reaction.gene_reaction_rule)
        gene_reaction_rule = or_re.sub("+", gene_reaction_rule)
        # To prevent shorter gene names from replacing substrings in
        # longer names, go in order from longest to shortest.
        reaction_genes = sorted(the_reaction._genes, reverse=True,
                                key=lambda x: len(x.id))
        # Replace each gene in the gpr string with 1 if it is still
        # active, or 0 if it is being knocked out.
        for gene in reaction_genes:
            if gene in gene_list:
                gene_reaction_rule = gene_reaction_rule.replace(gene.id, '0')
            else:
                gene_reaction_rule = gene_reaction_rule.replace(gene.id, '1')
        gene_reaction_rule = spontaneous_re.sub('1', gene_reaction_rule)
        if not eval(gene_reaction_rule):  # evaluates to 0 when gpr is false
            knocked_out_reactions.append(the_reaction)
    return knocked_out_reactions


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
