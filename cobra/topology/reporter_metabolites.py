# Based on Patil et al 2005 PNAS 102:2685-9
# TODO: Validate cobra.core compliance
from __future__ import print_function
from numpy import array, mean, std, where
from scipy.stats import norm, randint
from six import iteritems


def identify_reporter_metabolites(cobra_model, reaction_scores_dict,
                                  number_of_randomizations=1000,
                                  scoring_metric='default', score_type='p',
                                  entire_network=False,
                                  background_correction=True,
                                  ignore_external_boundary_reactions=False):
    """Calculate the aggregate Z-score for the metabolites in the model.
    Ignore reactions that are solely spontaneous or orphan. Allow the scores to
    have multiple columns / experiments.   This will change the way the output
    is represented.

    cobra_model: A cobra.Model object

    TODO: CHANGE TO USING DICTIONARIES for the_reactions: the_scores

    reaction_scores_dict:  A dictionary where the keys are reactions in
    cobra_model.reactions and the values are the scores.  Currently, only
    supports a single numeric value as the value; however, this will be updated
    to allow for lists

    number_of_randomizations: Integer.  Number of random shuffles of the
    scores to assess which are significant.

    scoring_metric: default means divide by k**0.5

    score_type: 'p' Is the only option at the moment and indicates p-value.

    entire_network: Boolean. Currently, only compares scores calculated from
    the_reactions

    background_correction: Boolean.  If True apply background correction to the
    aggreagate Z-score

    ignore_external_boundary_reactions: Not yet implemented. Boolean.  If True
    do not count exchange reactions when calculating the score.
    """

    # Add in a function to calculate based on correlation coefficients and to
    # deal with other multidimensional data.
    the_reactions = reaction_scores_dict.keys()
    the_scores = reaction_scores_dict.values()
    if score_type == 'p' and not hasattr(the_scores[0], '__iter__'):
        # minimum and maximum p-values are used to prevent numerical problems.
        # haven't decided whether an arbitrary min / max 1e-15 is preferred to
        # blunting the ends based on the values closest to 0 or 1.
        the_reactions = reaction_scores_dict.keys()
        the_scores = array(reaction_scores_dict.values())
        minimum_p = min(the_scores[the_scores.nonzero()[0]])
        maximum_p = max(the_scores[where(the_scores < 1)[0]])
        the_scores[where(the_scores < minimum_p)] = minimum_p
        the_scores[where(the_scores > maximum_p)] = maximum_p
        the_scores = -norm.ppf(the_scores)
        # update the dictionary with the new scores
        reaction_scores_dict = dict(zip(the_reactions, the_scores))
    elif hasattr(the_scores[0], '__iter__'):
        # In the case that the_scores is a list of lists, assume that each list
        # is the score for each reaction in the_reactions across all reactions.
        # Then for each metabolite, calculate the invnorm(|Pearson Correlation
        # Coefficient| for each reaction pair that it links.
        raise Exception("This isn't implemented yet")

    # Get the connectivity for each metabolite
    the_metabolites = set()
    for x in reaction_scores_dict:
        the_metabolites.update(x._metabolites)

    metabolite_scores = {}
    metabolite_connections = {}
    # Calculate the score for each metabolite
    for the_metabolite in the_metabolites:
        nonspontaneous_connections = [x for x in the_metabolite._reaction
                                      if x.gene_reaction_rule.lower() not in
                                      ['s0001', '']]
        tmp_score = 0
        number_of_connections = len(nonspontaneous_connections)
        for the_reaction in nonspontaneous_connections:
            if the_reaction not in reaction_scores_dict:
                if not entire_network:
                    number_of_connections -= 1
                continue
            else:
                tmp_score += reaction_scores_dict[the_reaction]
        metabolite_scores[the_metabolite] = tmp_score
        metabolite_connections[the_metabolite] = number_of_connections

    # NOTE: Doing the corrections based only on the significantly perturbed
    # scores is probably going to underestimate the significance.
    if background_correction:
        correction_dict = {}
        for i in set(metabolite_connections.values()):
            # if entire_network # add in a section to deal with the situation
            # where the entire network structure is considered by only have
            # p-values for a limited subset.
            #
            # Basically, what we're doing here is that for each i we select i
            # scores number_of_randomizations times
            the_random_indices = randint.rvs(
                0, len(the_scores), size=(number_of_randomizations, i))
            random_score_distribution = array(
                [sum(the_scores[x])
                 for x in list(the_random_indices)]) / i**0.5
            correction_dict[i] = [mean(random_score_distribution),
                                  std(random_score_distribution, ddof=1)]

    for the_metabolite, the_score in iteritems(metabolite_scores):
        number_of_connections = metabolite_connections[the_metabolite]
        if number_of_connections > 0:
            # Correct based on background distribution
            if background_correction:
                # if the list of scores is only for significant perturbations
                # then the background correction shouldn't be applied because
                # the current sampling method only takes into account
                # the_scores not the entire network.  It'd be more accurate to
                # assign unscored reactions a default score.
                the_score = ((the_score / number_of_connections**.5) -
                             correction_dict[number_of_connections][0]) / \
                    correction_dict[number_of_connections][1]
            else:
                the_score = the_score / number_of_connections**.5
            # Update the score
            metabolite_scores[the_metabolite] = the_score

    return_dictionary = {'scores': metabolite_scores,
                         'connections': metabolite_connections}
    if background_correction:
        return_dictionary['corrections'] = correction_dict

    return return_dictionary
