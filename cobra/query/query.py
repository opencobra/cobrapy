#cobra.query.query.py
#Will serve as a location to house the growing number of
#simple query functions attached to cobra.Model

#NOTE: Many of the find functions are gone because Reactions,
#Metabolites, and Genes are now away of each other.

import re
#####
def print_reactions_involving_metabolite(cobra_model, the_metabolites):
    """Update to allow for multiple metabolite search

    cobra_model: A cobra.Model object

    the_metabolites: A list of cobra.Metabolites or metabolite ids that are in
    cobra_metabolites.

    #TODO: Move this to the Metabolite class

    """
    if hasattr(the_metabolites, 'id'):
        the_metabolites = [the_metabolites]
    elif not hasattr(the_metabolites, '__iter__'):
        the_metabolites = [the_metabolites]
    if not hasattr(the_metabolites[0], 'id'):
        the_metabolites = [cobra_model.metabolites[cobra_model.metabolites.index(x)]
                           for x in the_metabolites]
        
    for the_metabolite in the_metabolties:
        for the_reaction in the_metabolite._reaction:
            print the_reaction.reaction
 
         
def get_translation_reactions(cobra_model, genes_of_interest):
    """Find the translation elongation reactions for a set of genes
    in a cobra model.  Related to ME-model extensions

    cobra_model:  A cobra.Model object.

    genes_of_interest:  A list of genes from cobra_model.genes.
    
    """
    gene_translation_reactions = defaultdict(list)
    for the_reaction in cobra_model.reactions:
        if 'translation_elongation' in the_reaction:
            for the_gene in genes_of_interest:
                if the_gene in the_reaction:
                    gene_translation_reactions[the_gene].append(the_reaction)
                    continue
    return gene_translation_reactions


if __name__ == '__main__':
    from cPickle import load
    from time import time
    solver = 'glpk'
    test_directory = '../test/data/'
    with open(test_directory + 'salmonella.pickle') as in_file:
        cobra_model = load(in_file)

    #TODO: Add in tests for each function
    print 'Need to add in tests for %s'%repr(['print_reactions_involving_metabolite'])
                                              
