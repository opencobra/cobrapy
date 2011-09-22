#cobra.manipulation.modify.py
import sys
#Here we're dealing with the fact that numpy isn't supported in jythong
#and we've made some numjy modules that use the cern.colt matrices to
#mimic the used numpy behavior.
if hasattr(sys, 'JYTHON_JAR'):
    raise Exception("Experimental modules of numpy/scipy for java that are" +\
    "not yet ready for prime time.")
    from cobra.scipy import sparse
    from cobra.numpy import where, array, ones
else:
    from scipy import sparse
    from numpy import where, array, ones
import re
from copy import deepcopy
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
        print 'All metabolites used in at least 1 reaction'

    
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
            print '%s not in %s'%(the_reaction.id,
                                      cobra_model.id)
    if not pruned_reactions:
        print 'All reactions have at least 1 metabolite'
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
    [setattr(x, 'functional', False)
     for x in gene_list]
    the_reactions = set()
    [the_reactions.update(x._reaction)
     for x in gene_list]
    for the_reaction in the_reactions:
        old_lower_bound = the_reaction.lower_bound
        old_upper_bound = the_reaction.upper_bound
        the_gene_reaction_relation = deepcopy(the_reaction.gene_reaction_rule)
        for the_gene in the_reaction._genes:
            the_gene_re = re.compile('(^|(?<=( |\()))%s(?=( |\)|$))'%re.escape(the_gene.id))
            if the_gene in gene_list:
                the_gene_reaction_relation = the_gene_re.sub('False', the_gene_reaction_relation)
            else:
                the_gene_reaction_relation = the_gene_re.sub('True', the_gene_reaction_relation)
        the_gene_reaction_relation = spontaneous_re.sub('True', the_gene_reaction_relation)
        if not eval(the_gene_reaction_relation):
            cobra_model._trimmed_reactions[the_reaction] = (old_lower_bound,
                                                            old_upper_bound)
            the_reaction.lower_bound = 0.
            the_reaction.upper_bound = 0.
            cobra_model._trimmed = True

    cobra_model._trimmed_genes =  list(set(cobra_model._trimmed_genes + gene_list))


def old_delete_model_genes(cobra_model, gene_list, gene_prefix='(STM|PSLT){1}',
                       cumulative_deletions=False, disable_orphans=False):
    """delete_model_genes will set the upper and lower bounds for reactions
    catalyzed by the genes in gene_list if deleting the genes means that
    the reaction cannot proceed according to
    cobra_model.reactions[:].gene_reaction_rule

    DEPRECATED

    TODO: Update this to use boolean logic based on the cobra_model.genes
    and cobra.Gene.  Once this is implemented, it will be much easier to
    perform this function and clean up the code.  It will also no longer
    be necessary to directly modify the bounds.

    gene_prefix: NO LONGER USED
    
    cumulative_deletions: False or True.  If True then any previous
    deletions will be maintained in the model.

    TODO: Rewrite to use dicts for _trimmed*

    TODO: All this will be replaced by Boolean logic associated with
    #the cobra.Gene.functional in cobra.Reaction.gene_reaction_rule

    TODO: Update this to refer to cobra.(Gene|Reaction) in the
    _trimmed_(genes|reactions) fields and remove _trimmed_indices
    
    """
    raise Exception('This does not work')
    if not hasattr(cobra_model, '_trimmed') or not cumulative_deletions:
        cobra_model._trimmed = False
        cobra_model._trimmed_genes = []
        cobra_model._trimmed_reactions = {} #Store the old bounds in here.
    spontaneous_re = re.compile('(^|(?<=( |\()))s0001(?=( |\)|$))')
    #Allow a single gene to be fed in as a string instead of a list.
    gene_list = deepcopy(gene_list)
    if not hasattr(gene_list, '__iter__') or \
           hasattr(gene_list, 'id'):  #cobra.Gene has __iter__
        gene_list = [gene_list]
    if hasattr(gene_list[0], 'id'):
        gene_list = [x.id for x in gene_list]
    #Allow gene_names or genes to be used
    if not gene_list[0] in cobra_model.genes:
        gene_names_in_model = [x.name for x in cobra_model.genes]
        gene_names = False
        tmp_list = []
        for the_gene in gene_list:
            #Deal with the case where a gene name maps to multiple loci
            the_count = gene_names_in_model.count(the_gene)
            tmp_index = 0
            for i in range(the_count):
                tmp_index = gene_names_in_model.index(the_gene, tmp_index)
                tmp_list.append(cobra_model.genes[tmp_index].id)
                gene_names = True
                tmp_index += 1
            if the_count > 1:
                print '%s in cobra_model.genes[:].name %i times'%(the_gene, the_count)
            elif the_count == 0:
                print '%s not in cobra_model.genes[:].name'%the_gene
        if not gene_names:
            print 'Cannot find %s in cobra_model.genes or cobra_model.genes[:].name'%repr(gene_list)
        gene_list = tmp_list
    #This deals with the fact that some genes have the same base but an
    #additional suffix
    #TODO: All this will be replaced by Boolean logic associated with
    #the cobra.Gene.functional in cobra.Reaction.gene_reaction_rule
    tmp_gene_reaction_dict = dict([(x, deepcopy(x.gene_reaction_rule))
                                       for x in cobra_model.reactions])
    reaction_to_index = dict(zip([x for x in cobra_model.reactions],
                                 range(len(cobra_model.reactions))))
    #Replace each gene from gene_list in the Boolean gene_protein_reaction
    #representations with False
    for the_gene in gene_list:
        the_gene_re = re.compile('(^|(?<=( |\()))%s(?=( |\)|$))'%re.escape(the_gene))
        for the_key, tmp_gpr in tmp_gene_reaction_dict.items():
            while the_gene_re.search(tmp_gpr):
                tmp_gpr = the_gene_re.sub('False', tmp_gpr)
            tmp_gene_reaction_dict[the_key] = tmp_gpr
    #All spontaneous items are changed to True
    tmp_gene_reaction_dict = dict([(k, spontaneous_re.sub('True', v))
                                   for k, v in tmp_gene_reaction_dict.items()])
    #Replace all genes not in gene_list with True.  Gene id's can
    #contain letters, digits, _, . and $.
    #There's probably a better regexp to do this but it hurts my head
    other_gene_re = re.compile('(?P<name>[a-zA-Z0-9_.$]+)')
    do_not_change_re = re.compile('^(( *\( *)+|(False)|(True)|(and)|(or)|| +|( *\) *)+)$')
    for the_key, tmp_gpr in tmp_gene_reaction_dict.items():
        new_gpr = ''
        for the_atom in other_gene_re.split(tmp_gpr):
            if do_not_change_re.match(the_atom):
                new_gpr += the_atom
            else:
                new_gpr += 'True'
        if new_gpr == '':
            if disable_orphans:
                new_gpr = 'False'
            else:
                new_gpr = 'True'
        tmp_gene_reaction_dict[the_key] = eval(new_gpr)

    #This goes based on reactions as being in a list.
    #Update to deal with cobra.Reaction objects
    cobra_model._trimmed_genes =  list(set(cobra_model._trimmed_genes + gene_list))
    for the_reaction, gpr_active in tmp_gene_reaction_dict.items():
        if not gpr_active:
            old_lower_bound = deepcopy(the_reaction.lower_bound)
            old_upper_bound = deepcopy(the_reaction.upper_bound)
            the_reaction.lower_bound = 0.
            the_reaction.upper_bound = 0.
            cobra_model._trimmed_reactions.update({the_reaction:
                                                   (old_lower_bound,
                                                    old_upper_bound)})
            cobra_model._trimmed = True

if __name__ == '__main__':
    from cPickle import load
    from time import time
    solver = 'glpk'
    test_directory = '../test/data/'
    with open(test_directory + 'salmonella.pickle') as in_file:
        cobra_model = load(in_file)

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
        print 'Successful deletion of %s'%repr(gene_list)
    else:
        print 'Failed deletion of %s\n%s reactions did not match'%(repr(gene_list),
                                                                   repr(symmetric_difference))






