#cobra.manipulation.modify.py
from copy import deepcopy
from .. import Reaction
from warnings import warn

def initialize_growth_medium(cobra_model, the_medium='MgM', 
                             external_boundary_compartment='e',
                             external_boundary_reactions=None,
                             reaction_lower_bound=0., reaction_upper_bound=1000.,
                             irreversible=False,
                             reactions_to_disable=None):
    """Sets all of the input fluxes to the model to zero and then will
    initialize the input fluxes to the values specified in the_medium if
    it is a dict or will see if the model has a composition dict and use
    that to do the initialization.

    cobra_model: A cobra.Model object.


    the_medium: A string, or a dictionary.  If a string then the
    initialize_growth_medium function expects that the_model has an
    attribute dictionary called media_compositions, which is a dictionary of
    dictionaries for various medium compositions.  Where a medium
    composition is a dictionary of external boundary reaction ids for the medium
    components and the external boundary fluxes for each medium component.


    external_boundary_compartment:  None or a string.  If not None then it specifies
    the compartment in which to disable all of the external systems boundaries.

    external_boundary_reactions: None or a list of external_boundaries that are to
    have their bounds reset.  This acts in conjunction with external_boundary_compartment.


    reaction_lower_bound: Float.  The default value to use for the lower
    bound for the boundary reactions.

    reaction_upper_bound: Float.  The default value to use for the upper
    bound for the boundary.

    irreversible: Boolean.  If the model is irreversible then the medium composition
    is taken as the upper bound

    reactions_to_disable: List of reactions for which the upper and lower bounds
    are disabled.  This is superceded by the contents of media_composition
  
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
    if external_boundary_reactions is not None:
        if isinstance(external_boundary_reactions[0], str):
            external_boundary_reactions = map(cobra_model.reactions.get_by_id, external_boundary_reactions)
    elif external_boundary_compartment is None:
            warn("We are initializing the medium without first adjusting all external boundary reactions")
        
    #Select the system_boundary reactions to reset
    if external_boundary_compartment is not None:
        _system_boundaries = dict([(x, x.get_compartments())
                                   for x in cobra_model.reactions
                                   if x.boundary == 'system_boundary'])
        [_system_boundaries.pop(k) for k, v in list(_system_boundaries.items())
         if len(v) == 1 and external_boundary_compartment not in v]
        if external_boundary_reactions is None:
            external_boundary_reactions = _system_boundaries.keys()
        else:
            external_boundary_reactions += _system_boundaries.keys()


    for the_reaction in external_boundary_reactions:
        the_reaction.lower_bound = reaction_lower_bound
        if the_reaction.upper_bound == 0:
            the_reaction.upper_bound = reaction_upper_bound
    #Disable specified reactions
    if reactions_to_disable is not None:
        if isinstance(reactions_to_disable[0], str):
            reactions_to_disable = map(cobra_model.reactions.get_by_id, reactions_to_disable)
        for the_reaction in reactions_to_disable:
            the_reaction.lower_bound = the_reaction.upper_bound = 0.


    #Update the model inputs based on the_medium
    for the_component in medium_composition.keys():
        the_reaction = cobra_model.reactions.get_by_id(the_component)
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
            reaction_dict = dict([(k, v*-1)
                                  for k, v in reaction._metabolites.items()])
            reverse_reaction.add_metabolites(reaction_dict)
            reverse_reaction._model = reaction._model
            reverse_reaction._genes = reaction._genes
            for gene in reaction._genes:
                gene._reaction.add(reverse_reaction)
            reverse_reaction._gene_reaction_rule = reaction._gene_reaction_rule
            reactions_to_add.append(reverse_reaction)
    cobra_model.add_reactions(reactions_to_add)
 
def revert_to_reversible(cobra_model, update_solution=True):
    """This function will convert a reversible model made by convert_to_irreversible
    into a reversible model.

    cobra_model:  A cobra.Model object which will be modified in place.

    NOTE: It might just be easiest to include this function in the Reaction class
    
    """
    reverse_reactions = [x for x in cobra_model.reactions
                         if x.reflection is not None and
                         x.id.endswith('_reverse')]

    for reverse in reverse_reactions:
        forward = reverse.reflection
        forward.lower_bound = -reverse.upper_bound
        forward.reflection = None
    #Since the metabolites and genes are all still in
    #use we can do this faster removal step.  We can
    #probably speed things up here.
    cobra_model.remove_reactions(reverse_reactions)
    # fix the solution
    if update_solution and cobra_model.solution is not None:
        x_dict = cobra_model.solution.x_dict
        for reverse in reverse_reactions:
            forward = reverse.reflection
            x_dict[forward.id] -= x_dict.pop(reverse.id)
        cobra_model.solution.x = [x_dict[r.id] for r in cobra_model.reactions]

   
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
