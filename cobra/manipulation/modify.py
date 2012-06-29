#cobra.manipulation.modify.py
from copy import deepcopy
from cobra.core.Reaction import Reaction
def add_exchange_reaction(cobra_model, the_metabolites,
                          reaction_reversibility=False):
    """Adds exchange reactions to a model for a set of metabolites.

    cobra_model: A Model object.

    the_metabolites: A cobra.Metabolite or list of cobra.Metabolites

    reaction_reversibility: True or False.  Indicates whether the reactions
    should be reversible.


    #TODO:  This is not compliant with current 
    """
    if not isinstance(the_metabolites, list):
        the_metabolites = [the_metabolites]
    if reaction_reversibility:
        lower_bound = -1000
    else:
        lower_bound = 0
    the_reactions = []
    for the_metabolite in the_metabolites:
        the_reaction = Reaction('EX_' + the_metabolite.id)
        the_reaction.name = the_metabolite + ' exchange reaction'
        the_reaction.add_metabolites({the_metabolite: -1})
        the_reaction.reversibility = reaction_reversibility
        the_reaction.lower_bound = lower_bound
        the_reaction.upper_bound = 1000
        the_reaction.objective_coefficient = 0.
        the_reaction.boundary = True
        the_reactions.append(the_reaction)
    cobra_model.add_reactions(the_reactions)
    

def initialize_growth_medium(cobra_model, the_medium='MgM', exchange_reactions=None,
                             reactions_to_disable=None,
                             exchange_lower_bound=0., exchange_upper_bound=1000.,
                             irreversible=False):
    """Sets all of the input fluxes to the model to zero and then will
    initialize the input fluxes to the values specified in the_medium if
    it is a dict or will see if the model has a composition dict and use
    that to do the initialization.

    cobra_model: A cobra.Model object.
    
    the_medium: Is a string, or a dictionary.  If a string then the
    initialize_growth_medium function expects that cobra_model has an
    attribute dictionary called media_compositions, which is a dictionary of
    dictionaries for various medium compositions.  Where a medium
    composition is a dictionary of exchange reaction ids for the medium
    components and the exchange fluxes for each medium component; note that
    these fluxes must be negative because they are being exchanged into the
    system

    exchange_reactions: None or a list of exchange reactions.  Because not all
    exchange reactions are required to start with EX_ in all models allow
    the user to specify the reactions

    exchange_lower_bound: Float.  The default value to use for the lower
    bound for the exchange reactions.

    exchange_upper_bound: Float.  The default value to use for the upper
    bound for the exchange_reactions.

    irreversible: Boolean.  If the model is irreversible then the medium composition
    is taken as the upper bound
  
    """
    #Zero all of the inputs to the model
    if hasattr(the_medium, 'keys'):
        medium_composition = the_medium
    else:
        if hasattr(cobra_model, 'media_compositions'):
            if the_medium in cobra_model.media_compositions:
                medium_composition = cobra_model.media_compositions[the_medium]
            else:
                raise Exception("%s is not in the model's media list"%the_medium)     
        else:
            raise Exception("the model doesn't have attribute media_compositions and the medium is not a dict")
    if exchange_reactions is None:
        #NOTE this is a bad way to do things.  Should just find reactions that cross
        #whatever the external boundary
        exchange_reactions = [x for x in cobra_model.reactions if x.startswith('EX_')]
    if not hasattr(exchange_reactions[0], 'id'):
        exchange_indices = map(cobra_model.reactions.index, exchange_reactions)
        #SPEED UP
        exchange_reactions = [cobra_model.reactions[x]
                              for x in exchange_indices]
    for the_reaction in exchange_reactions:
        the_reaction.lower_bound = exchange_lower_bound
        if the_reaction.upper_bound == 0:
            the_reaction.upper_bound = exchange_upper_bound
    if reactions_to_disable is not None:
        disabled_indices = map(cobra_model.reactions.index, reactions_to_disable)
        disabled_reactions = [cobra_model.reactions[x]
                              for x in disabled_indices]
        [(setattr(x, 'lower_bound', 0.),
          setattr(x, 'upper_bound', 0.))
         for x in disabled_reactions]

    #Update the model inputs based on the_medium
    for the_component in medium_composition.keys():
        the_index = cobra_model.reactions.index(the_component)
        the_reaction = cobra_model.reactions[the_index]
        if irreversible:
            the_reaction.upper_bound = medium_composition[the_component]
        else:
            the_reaction.lower_bound = medium_composition[the_component]


def convert_to_irreversible(cobra_model):
    """Will break all of the reversible reactions into two separate irreversible reactions with
    different directions.  Useful for some modeling problems.

    cobra_model: A Model object which will be modified in place.

    
    TODO: Can we just use a -1*guided_copy or something else?
    """
    reactions_to_add = []
    for reaction in cobra_model.reactions:
        #Potential bug because a reaction might run backwards naturally
        #and this would result in adding an empty reaction to the
        #model in addition to the reverse reaction.
        if reaction.lower_bound < 0:
            reverse_reaction = Reaction(reaction.id + "_reverse")
            reverse_reaction.lower_bound = 0
            reverse_reaction.upper_bound = reaction.lower_bound * -1
            reaction.lower_bound = 0
            #Make the directions aware of each other
            reaction.reflection = reverse_reaction
            reverse_reaction.reflection = reaction
            reaction.reversibility = reverse_reaction.reversibility = 0
            reaction_dict = dict([(k, v*-1)
                                  for k, v in reaction._metabolites.items()])
            reverse_reaction.add_metabolites(reaction_dict)
            reverse_reaction._model = reaction._model
            reverse_reaction._genes = reaction._genes
            reverse_reaction.gene_reaction_rule = reaction.gene_reaction_rule
            reactions_to_add.append(reverse_reaction)
    cobra_model.add_reactions(reactions_to_add)
 
def revert_to_reversible(cobra_model):
    """This function will convert a reversible model made by convert_to_irreversible
    into a reversible model.

    cobra_model:  A cobra.Model object which will be modified in place.

    NOTE: It might just be easiest to include this function in the Reaction class
    
    """
    reversible_reactions = [x for x in cobra_model.reactions
                            if x.reflection is not None and
                            not x.id.endswith('_reverse')]

    for the_reaction in reversible_reactions:
        the_reflection = the_reaction.reflection
        the_reaction.lower_bound = -the_reflection.lower_bound
        the_reaction.reflection = None
        #Since the metabolites and genes are all still in
        #use we can do this faster removal step.  We can
        #probably speed things up here.
        cobra_model.reactions.remove(the_reaction)

   
def convert_rule_to_boolean_rule(cobra_model, the_rule,
                                 return_gene_indices=False,
                                 index_offset=0):
    """Used to convert a text based gpr to an index based gpr.  This
    will also update the cobra_model.

    the_rule: A COBRA 2.0 compliant GPR string

    return_gene_indices: Boolean return the indices for the genes

    index_offset: Integer.  Set to 1 if the rules need to be generated for
    base 1 software like MATLAB.

    TODO: Test now that cobra.Gene is in use
    DEPRECATED:  This should be moved to the mlab module
    
    """
    raise Exception('cobra.manipulation.modify.convert_rule_to_boolean_rule is no longer functional')
    if the_rule == '':
        boolean_rule = ''
        the_gene_indices = []
    else:
        #POTENTIAL BUG
        #Deal with Simpheny style GPRs.  Can we use re or will it be too problematic?
        the_rule = the_rule.replace('+', ' and ').replace(',', ' or ')
        the_gene_list = the_rule.replace('(','').replace(')','').replace(' and ', '\t').replace(' or ','\t').replace(' ','').split('\t') 
        the_gene_indices = []
        new_rule = ''
        for the_gene in the_gene_list:
            #Add the gene to the model if it isn't there
            if not isinstance(the_gene, Gene):
                the_gene = Gene(the_gene)
            if the_gene.id in cobra_model.genes:
                the_gene_index = cobra_model.genes.index(the_gene.id) + index_offset
            else:
                the_gene_index = len(cobra_model.genes) + index_offset
                cobra_model.genes.append(the_gene)
            the_gene_indices.append(the_gene_index)
            #Split the_rule based on the first occurrence of the_gene
            #and set the_rule to the_suffix so that it can be split
            #for the next gene.  This will prevent problems with genes
            #being substrings of other genes or of gene indices.
            the_prefix, the_rule = the_rule.split(the_gene.id, 1)
            new_rule =  '%s%sx(%i)'%(new_rule, the_prefix, the_gene_index)
        the_rule = new_rule + the_rule
        boolean_rule = the_rule.replace(' and ', ' & ').replace(' or ', ' | ')
    if return_gene_indices:
        return {'boolean_rule':boolean_rule, 'gene_indices':the_gene_indices}
    else:
        return boolean_rule

if __name__ == '__main__':
    from time import time
    from cobra.test import create_test_model
    cobra_model = create_test_model()


    print 'Move this to test'
    irreversible_model = deepcopy(cobra_model)
    start_time = time()
    convert_to_irreversible(irreversible_model)
    print 'Convert to irreversible took %1.3f seconds'%(time()-start_time)

    reverted_model = deepcopy(irreversible_model)
    start_time = time()
    revert_to_reversible(reverted_model)
    print 'Convert to reversible took %1.3f seconds'%(time()-start_time)
