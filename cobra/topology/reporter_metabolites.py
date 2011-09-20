#cobra.topology.reporter_metabolites.py: Module for topological analysis of cobra_models
#Based on Patil et al 2005 PNAS 102:2685-9
#TODO: Validate cobra.core compliance
from copy import deepcopy
from numpy import array, corrcoef, mean, std, tril, where, unique, zeros
from scipy.stats import norm, randint
from collections import defaultdict

def identify_reporter_metabolites(cobra_model, the_reactions, the_scores,
                                  number_of_randomizations=1000, number_of_layers=1,
                                  scoring_metric='default', score_type='p',
                                  entire_network=False, background_correction=True,
                                  return_all=False, ignore_exchange_reactions=True):
    """Calculate the aggregate Z-score for the metabolites in the model.
    Ignore reactions that are solely spontaneous or orphan. Allow the scores to
    have multiple columns / experiments.   This will change the way the output
    is represented.

    cobra_model: A cobra.Model object

    the_reactions: A list of reactions for which to calculate aggregate scores.

    the_scores: Corresponding scores for the_reactions of type score_type.

    number_of_randomizations: Integer.  Number of random shuffles of the
    scores to assess which are significant.

    number_of_layers: 1 is the only option supported
    
    scoring_metric:  default means divide by k**0.5

    score_type: 'p' Is the only option at the moment and indicates p-value.

    entire_network: Boolean.  Currently, only compares scores calculated from the_reactions

    background_correction: Boolean.  If True apply background correction to the aggreagate
    Z-score

    return_all: Boolean.  If True return all scores.  If False only return scores > 0.

    ignore_exchange_reactions: Boolean.  If True do not count exchange reactions when
    calculating the score.

    #TODO: Rewrite to take advantage of cobra.Objects

    
    """
    cobra_model.update_stoichiometric_matrix()
    #Add in a function to calculate based on correlation coefficients and to
    #deal with other multidimensional data. 
    if score_type == 'p' and not hasattr(the_scores[0], '__iter__'):
        #minimum and maximum p-values are used to prevent numerical problems.
        #haven't decided whether an arbitrary min / max 1e-15 is preferred to
        #blunting the ends based on the values closest to 0 or 1.
        the_scores = array(the_scores)
        minimum_p = min(the_scores[ the_scores.nonzero()[0]])
        maximum_p = max(the_scores[where(the_scores < 1)[0]])
        the_scores[where(the_scores < minimum_p)] = minimum_p
        the_scores[where(the_scores > maximum_p)] = maximum_p
        the_scores = -norm.ppf(the_scores)
    elif hasattr(the_scores[0], '__iter__'):
        #In the case that the_scores is a list of lists, assume that each list is
        #the score for each reaction in the_reactions across all reactions.  Then
        #for each metabolite, calculate the invnorm(|Pearson Correlation
        #Coefficient| for each reaction pair that it links.
        tmp_scores = []
        for the_row in cobra_model._S.rows:
            tmp_reaction_indices = [the_reactions.index(y)
                                     for y in [cobra_model.reactions[x]
                                                for x in the_row] if y in the_reactions ]
            if len(tmp_reaction_indices) > 1:
                tmp_array = tril(abs(corrcoef([the_scores[x]
                                                  for x in tmp_reaction_indices])), -1)
                #HERE 2010-03-16
                #TODO: Make sure this is correct and then deal with the latter
                #part of the function calculating for p-values
                tmp_p_sum = tmp_array.sum()
                tmp_k = (tmp_array.shape[0]*(tmp_array.shape[0]-1))/2
                tmp_score = tmp_p_sum / tmp_k**.5
                tmp_scores.append(tmp_score)
            else:
                tmp_scores.append(0)
        the_scores = tmp_scores
    #Get the indices of spontaneous reactions.  They should not be counted as connections.
    nonenzymatic_reaction_indices = [cobra_model.reactions.index(x)
                                     for x in cobra_model.reactions 
                                     if x.gene_reaction_rule == 's0001' or
                                     x.gene_reaction_rule == '']
    if hasattr(the_reactions[0], 'id'):
        scoring_dict = dict(zip(map(lambda x: x.id, the_reactions), the_scores))
    else:
        scoring_dict = dict(zip(the_reactions, the_scores))
    #Create a matrix of the stoichiometry for only the reactions involved in
    #the reactions.
    reaction_indices = array(map(cobra_model.reactions.index, the_reactions))
    #A fast way to identify the connectivity distribution for the metabolites
    #need to have a function to feed to the defulat dict.
    def zero():
        return(0)
    connectivity_dict = defaultdict(zero)
    #Here, we transpose the S matrix because it's faster to access rows in
    #the scipy sparse matrix.  Then we use the rows function which will provide
    #the indice of each column (metabolite now that the matrix has been transposed)
    #that is nonzero.  Then we select out the rows (reactions) that correspond to
    #a reaction in the reactions.
    #To get the connectivity for each metabolite, we just count the number of
    #times the index is detected in a row.
    for x in reduce(lambda x, y: x+y,cobra_model._S.T.rows[reaction_indices]):
        connectivity_dict[x] += 1
    #Here we just get the distinct connectivity values.
    connectivity_values = list(set(connectivity_dict.values()))

    correction_dict = {}
    #NOTE: Doing the corrections based only on the significantly perturbed scores
    #is probably going to underestimate the significance.
    if background_correction:
        for i in connectivity_values:
            #if entire_network # add in a section to deal with the situation where
            #the entire network structure is considered by only have p-values for
            #a limited subset.
            #
            #Basically, what we're doing here is that for each i we select i
            #scores number_of_randomizations times
            the_random_indices = randint.rvs(0,len(the_scores),size=(number_of_randomizations, i))

            random_score_distribution = array([sum(the_scores[x]) for x in list(the_random_indices)]) /i**0.5
            correction_dict[i] = [mean(random_score_distribution),
                                      std(random_score_distribution)] 

    metabolite_scores = []
    metabolite_connections = []
    #Calculate the score for each metabolite
    for the_row in cobra_model._S.rows:
        #Currently, we do all the metabolites in model, but we could speed this
        #up by only dealing with the keys from connectivity_dict
        k = 0
        tmp_score = 0
        for reaction_index in the_row:
            the_reaction = cobra_model.reactions[reaction_index]
            if hasattr(the_reaction, 'id'):
                the_reaction = the_reaction.id
            if the_reaction in scoring_dict.keys():
                tmp_score += scoring_dict[the_reaction]
                k += 1
            elif entire_network and reaction_index not in nonenzymatic_reaction_indices:
                #This increases the k for reactions that are not in the_reactions
                #but exist in cobra_model and have an associated enzyme.
                k += 1
        if k > 0:
            #Correct based on background distribution
            if background_correction:
                #if the list of scores is only for significant perturbations then the
                #background correction shouldn't be applied because the current sampling
                #method only takes into account the_scores not the entire network.
                #It'd be more accurate to assign unscored reactions a default score.
                tmp_score = ((tmp_score / k**.5)-correction_dict[k][0])/correction_dict[k][1]
            else:
                tmp_score = tmp_score / k**.5
        metabolite_scores.append(tmp_score)
        metabolite_connections.append(k)
    if entire_network or return_all:
        if background_correction:
            return([deepcopy(cobra_model.metabolites),  array(metabolite_scores),
                    array(metabolite_connections) , correction_dict])
        else:
            return([deepcopy(cobra_model.metabolites), array(metabolite_scores),
                    array(metabolite_connections)])
    else:
        the_indices = list(where(array(metabolite_scores) > 0)[0])
        return([[cobra_model.metabolites[x] for x in the_indices],
                array([metabolite_scores[x] for x in the_indices]),
                array([metabolite_connections[x] for x in the_indices])])

def ppmap_identify_reporter_metabolites(keywords):
    """
    A function that receives a dict with all of the parameters for identify_reporter_metabolites
    Serves to make it possible to call the reporter metabolites function from ppmap.
    It only will be useful for parallel experiments not for breaking up a single experiment.
    
    """
    the_results = identify_reporter_metabolites(**keywords)
    return({'id': the_id, 'results': the_results })

if __name__ == '__main__':
    from cPickle import load
    from time import time
    solver = 'glpk'
    test_directory = '../test/data/'
    with open(test_directory + 'salmonella.pickle') as in_file:
        cobra_model = load(in_file)
    with open(test_directory + 'reaction_p_values.pickle') as in_file:
        reaction_p = load(in_file)

    the_reactions = reaction_p.keys()
    the_scores = [reaction_p[k] for k in the_reactions]
    tmp_reps = identify_reporter_metabolites(cobra_model, the_reactions,
                                             the_scores,
                                             background_correction=True,
                                             return_all=True)

    print 'Need to add in validation for the test'
