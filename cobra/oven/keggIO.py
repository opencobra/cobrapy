import csv
import re
import copy
from os.path import join, abspath, dirname

import cobra

# the default file locations
kegg_directory = join(abspath(dirname(__file__)), "kegg_files")
keggdictpath_default = join(kegg_directory, "kegg_dict.csv")
reactionlst_default = join(kegg_directory, "reaction.lst")
blacklistpath_default = join(kegg_directory, "kegg_blacklist.csv")


def intify(string):
    """returns integer representation of the str
    If str is a single letter, it will return 1"""
    if string.isdigit():
        return int(string)
    # if the expression contains n, the default value is 2
    n = 2
    if string == "2n":
        return 2 * n
    try:
        return eval(string)
    except:
        raise ValueError(string)


def parse_split_array(str_array):
    """takes in an array of strings, each of which is either
    - a compound OR
    - a number followed by a compound
    returns [array_of_metabolites, corresponding_coefficient]"""
    metabolites = []
    coefficients = []
    for string in str_array:
        string = string.strip()
        if string[0].isupper():  # starts with an uppercase letter
            # there is no number associated, so it should be 1
            metabolites.append(string)
            coefficients.append(1)
        else:
            the_coefficient, the_metabolite = string.split()
            metabolites.append(the_metabolite)
            coefficients.append(intify(the_coefficient))
    return [metabolites, coefficients]


def import_kegg_reactions(compartment="c", reactionlstpath=None,
                        keggdictpath=None, blacklistpath=None):
    """reads in kegg reactions from the three given files
    compartment: the compartment to which each reaction will be added

    If no file is specified for any of these, a default file will be used:
    reactionlstpath: path to path of kegg reactions
        the format should be
        reactionid: Met1 + 2 Met2 <=> Met3 + 2 Met4
    keggdictpath: path to csv file translating between kegg and cobra
        metabolite ID's
        first column contains kegg ID, second contains cobra id
    blacklistpath: path to csv containing the kegg blacklist
        first columm contains the kegg id's of blacklisted reactions

    returns: cobra model with all of the included reactions"""

    if reactionlstpath is None:
        reactionlstpath = reactionlst_default
    if keggdictpath is None:
        keggdictpath = keggdictpath_default
    if blacklistpath is None:
        blacklistpath = blacklistpath_default

    # read in kegg dictionary to translate between kegg and cobra id's
    keggdictfile = open(keggdictpath, "r")
    keggdictcsv = csv.reader(keggdictfile)
    keggdict = {}
    for line in keggdictcsv:
        keggdict[line[1]] = line[0]
    keggdictfile.close()
    # read in the kegg blacklist
    keggblacklistfile = open(blacklistpath, "r")
    keggblacklistcsv = csv.reader(keggblacklistfile)
    keggblacklist = []
    for line in keggblacklistcsv:
        keggblacklist.append(line[0])
    keggblacklistfile.close()

    # parse the file of kegg reactions
    keggfile = open(reactionlstpath, "r")
    # regular expressions to split strings
    colon_sep = re.compile(":").split
    arrow_sep = re.compile("<=>").split
    plus_sep = re.compile(" \+ ").split
    keggreactions = []
    cobra_reactions = []
    used_metabolites = {}
    for line in keggfile:
        [id, reactionstr] = colon_sep(line, maxsplit=1)
        # remove whitespace
        id = id.strip()
        # if the id is in the blacklist, no need to proceed
        if id in keggblacklist:
            continue
        # split into reactants and products
        reactants_str, products_str = arrow_sep(reactionstr, maxsplit=1)
        # break up reactant and product strings into arrays of
        # metabolites and coefficients
        reactant_metabolites, reactant_coefficients = \
            parse_split_array(plus_sep(reactants_str))
        product_metabolites, product_coefficients = \
            parse_split_array(plus_sep(products_str))
        # reactant coefficients all need to be multiplied by -1
        for i, coeff in enumerate(reactant_coefficients):
            reactant_coefficients[i] = coeff * -1
        # make one array for all compoenents
        kegg_metabolites = reactant_metabolites
        coefficients = reactant_coefficients
        kegg_metabolites.extend(product_metabolites)
        coefficients.extend(product_coefficients)
        # translate the metabolites from kegg to cobra
        metabolites = []
        try:
            for the_kegg_metabolite in kegg_metabolites:
                metabolites.append(keggdict[the_kegg_metabolite])
        # if one of the metabolites is not found, skip to the next line
        except KeyError:
            continue

        # make a Kegg reaction
        reaction = cobra.Reaction(id)
        metabolite_dict = {}  # dict of {metabolite : coefficient}
        for i, the_metabolite in enumerate(metabolites):
            metabolite_id = the_metabolite + "_" + compartment
            # if the metabolite already exists
            if metabolite_id in used_metabolites:
                used_metabolites[metabolite_id] = coefficients[i]
            else:
                # use a new metabolite
                new_metabolite = cobra.Metabolite(metabolite_id)
                used_metabolites[metabolite_id] = new_metabolite
                metabolite_dict[cobra.Metabolite(metabolite_id)] = \
                    coefficients[i]
        reaction.add_metabolites(metabolite_dict)
        reaction.notes["temporary_gapfilling_type"] = "Universal"
        # because the model will be converted to irreversible
        reaction.lower_bound = -1 * reaction.upper_bound
        cobra_reactions.append(reaction)
    keggfile.close()
    # add all of the reactions to a cobra model
    Universal = cobra.Model("Kegg_Universal_Reactions")
    Universal.add_reactions(cobra_reactions)
    return Universal
if __name__ == "__main__":
    from time import time
    start_time = time()
    test_import = import_kegg_reactions()
    duration = time() - start_time
    print "imported %d reactions in %.2f sec" % \
        (len(test_import.reactions), duration)
