#cobra.manipulation.modify.py
import re
from copy import deepcopy

from warnings import warn

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

    reactions_to_prune: None, a string matching a reaction.id, a cobra.Reaction, or
    as list of the ids / Reactions to remove from cobra_model.
    If None then the function will delete reactions that
    have no active metabolites in the model.

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
    """Undoes the effects of a call to delete_model_genes. Modifies cobra_model in place.

    cobra_model:  A cobra.Model object

    """

    if cobra_model._trimmed_genes is not None:
        [setattr(x, 'functional', True)
         for x in cobra_model._trimmed_genes]
    
    if cobra_model._trimmed_reactions is not None:
        for the_reaction, (lower_bound,
                           upper_bound) in cobra_model._trimmed_reactions.items():
            the_reaction.lower_bound = lower_bound
            the_reaction.upper_bound = upper_bound
    #
    cobra_model._trimmed_genes = None
    cobra_model._trimmed_reactions = None
    cobra_model._trimmed = False
    #Reset these to deal with potential bugs from users accessing private variables
    for the_attribute in  ['_lower_bounds', '_upper_bounds',
                           '_S', '_objective_coefficients']:
        if hasattr(cobra_model, the_attribute):
            setattr(cobra_model, the_attribute, None)

spontaneous_re = re.compile('(^|(?<=( |\()))s0001(?=( |\)|$))')
def find_gene_knockout_reactions(cobra_model, gene_list):
    """identify reactions which will be disabled when the genes are knocked out"""

    potential_reactions = set()
    for x in gene_list:
        potential_reactions.update(x._reaction)

    knocked_out_reactions = []
    for the_reaction in potential_reactions:
        # operate on a copy
        gene_reaction_rule = "" + the_reaction.gene_reaction_rule
        for the_gene in the_reaction._genes:
            if the_gene in gene_list:
                gene_reaction_rule = gene_reaction_rule.replace(the_gene.id, 'False')
            else:
                gene_reaction_rule = gene_reaction_rule.replace(the_gene.id, 'True')
        gene_reaction_rule = spontaneous_re.sub('True', gene_reaction_rule)
        if not eval(gene_reaction_rule):
            knocked_out_reactions.append(the_reaction)
    return knocked_out_reactions


def delete_model_genes(cobra_model, gene_list,
                       cumulative_deletions=False, disable_orphans=False):
    """delete_model_genes will set the upper and lower bounds for reactions
    catalyzed by the genes in gene_list if deleting the genes means that
    the reaction cannot proceed according to
    cobra_model.reactions[:].gene_reaction_rule

    cumulative_deletions: False or True.  If True then any previous
    deletions will be maintained in the model.

    TODO: Rewrite to use dicts for _trimmed*

    TODO: All this will be replaced by Boolean logic associated with
    #the cobra.Gene.functional in cobra.Reaction.gene_reaction_rule

    TODO: Update this to refer to cobra.(Gene|Reaction) in the
    _trimmed_(genes|reactions) fields and remove _trimmed_indices
    
    """
    if not hasattr(cobra_model, '_trimmed') or not cumulative_deletions:
        cobra_model._trimmed = False
        cobra_model._trimmed_genes = []
        cobra_model._trimmed_reactions = {} #Store the old bounds in here.
    spontaneous_re = re.compile('(^|(?<=( |\()))s0001(?=( |\)|$))')
    #Allow a single gene to be fed in as a string instead of a list.
    if not hasattr(gene_list, '__iter__') or \
           hasattr(gene_list, 'id'):  #cobra.Gene has __iter__
        gene_list = [gene_list]

    if not hasattr(gene_list[0], 'id'):
        if gene_list[0] in cobra_model.genes:
                tmp_gene_dict = dict([(x.id, x) for x in cobra_model.genes])
        else:
            #assume we're dealing with names if no match to an id
            tmp_gene_dict = dict([(x.name, x) for x in cobra_model.genes])
        gene_list = [tmp_gene_dict[x] for x in gene_list]

    #Make the genes non-functional
    for x in gene_list:
        x.functional = False

    for the_reaction in find_gene_knockout_reactions(cobra_model, gene_list):
        old_lower_bound = the_reaction.lower_bound
        old_upper_bound = the_reaction.upper_bound
        cobra_model._trimmed_reactions[the_reaction] = (old_lower_bound,
                                                        old_upper_bound)
        the_reaction.lower_bound = 0.
        the_reaction.upper_bound = 0.
        cobra_model._trimmed = True

    cobra_model._trimmed_genes =  list(set(cobra_model._trimmed_genes + gene_list))



if __name__ == '__main__':
    from time import time
    from cobra.test import create_test_model
    cobra_model = create_test_model()


    #TODO: Add in tests for each function
    cumulative_deletions=False
    disable_orphans=False
    gene_list = ['STM1067', 'STM0227']
    #The following reactions are trimmed when STM1332 and STM1101 are deleted
    dependent_reactions = set(['3HAD121',
                               '3HAD160',
                               '3HAD80',
                               '3HAD140',
                               '3HAD180',
                               '3HAD100',
                               '3HAD181',
                               '3HAD120',
                               '3HAD60',
                               '3HAD141',
                               '3HAD161',
                               'T2DECAI',
                               '3HAD40'])
    delete_model_genes(cobra_model, gene_list)
    symmetric_difference = dependent_reactions.symmetric_difference([x.id for x in cobra_model._trimmed_reactions])
    if len(symmetric_difference) == 0:
        'Successful deletion of %s'%repr(gene_list)
    else:
        'Failed deletion of %s\n%s reactions did not match'%(repr(gene_list),
                                                                   repr(symmetric_difference))






