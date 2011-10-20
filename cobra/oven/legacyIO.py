from os.path import isfile
from csv import reader
from re import compile, findall
from warnings import warn

from scipy.io import loadmat, savemat

import cobra
from cobra.io.sbml import create_cobra_model_from_sbml_file


def _fix_legacy_id(the_id):
    the_id = the_id.replace('_DASH_', '__')
    the_id = the_id.replace('_FSLASH_', '/')
    the_id = the_id.replace('_BSLASH_', "\\")
    the_id = the_id.replace('_LPAREN_', '(')
    the_id = the_id.replace('_LSQBKT_', '[')
    the_id = the_id.replace('_RSQBKT_', ']')
    the_id = the_id.replace('_RPAREN_', ')')
    the_id = the_id.replace('_COMMA_', ',')
    the_id = the_id.replace('_PERIOD_', '.')
    the_id = the_id.replace('_APOS_', "'")
    the_id = the_id.replace('&amp;', '&')
    the_id = the_id.replace('&lt;', '<')
    the_id = the_id.replace('&gt;', '>')
    the_id = the_id.replace('&quot;', '"')
    return the_id


def read_legacy_sbml(filename):
    model = create_cobra_model_from_sbml_file(filename)
    for reaction in model.reactions:
        reaction.id = _fix_legacy_id(reaction.id)
    for metabolite in model.metabolites:
        metabolite.id = _fix_legacy_id(metabolite.id)
    return model


def _header_count(filename):
    """count the number of header lines in a file
    The header is defined as over when a line is found which begins
    with a number"""
    file = open(filename, "r")
    for i, line in enumerate(file):
        if line[0].isdigit():
            file.close()
            return i
    file.close()
    return False


def _open_and_skip_header(filename):
    """returns a csv file with the header skipped"""
    count = _header_count(filename)
    if not count:
        raise (IOError, "%s corrupted" % filename)
    file = open(filename, "r")
    for i in range(count):
        file.readline()
    return reader(file, delimiter="\t")


def _find_metabolites_by_base(base, metabolites):
    """search for all metabolites in the list which match the base name.
    For example, "h2o" will identify both "h2o(c)" and "h2o(e)" """
    search = compile(base + "\([a-z]\)")
    found = []
    for the_metabolite in metabolites:
        if search.match(the_metabolite.id) is not None:
            found.append(the_metabolite)
    return found


def read_simpheny(baseName, minlowerbound=-1000, maxupperbound=1000,
        maximize_info=True):
    """Imports SimPheny Exported files. This is a legacy import, and
    may give an erroenous model. Use with caution.

    baseName: The filepath to the exported SimPheny files without any of
    the extensions. On Windows, it helps if baseName is a raw string
    (i.e. r"Path\\to\\files")

    maximize_info: An optional boolean keyword argument. If True, then an
    attempt will be made to parse the gpr and metabolite info files, and the
    function will take a little bit longer.
    NOTE: Simpheny gene protein relationships have been known to be corrupted

    returns: the imported simpheny model"""

    # check to make sure the files can be read
    if not(isfile(baseName + ".met")
        and isfile(baseName + ".rxn")
        and isfile(baseName + ".sto")):
        # try again with modifying the baseName
        baseName = baseName.encode("string-escape")
        if not(isfile(baseName + ".met")
            and isfile(baseName + ".rxn")
            and isfile(baseName + ".sto")):
            raise (IOError, "Input file(s) not found")
    model = cobra.Model("SimPheny import from " + baseName)
    # read in metabolite file
    metfile = _open_and_skip_header(baseName + ".met")
    metabolites = []
    for line in metfile:
        if len(line) == 0:
            break
        metabolite = cobra.Metabolite(id=line[1], name=line[2],
            compartment=line[3])
        if maximize_info:
            compartment_search = findall("\([a-z]\)$", metabolite.id)
            if compartment_search != []:
                metabolite.compartment = compartment_search[0][1]
                model.compartments[metabolite.compartment] = line[3]
        metabolites.append(metabolite)
    model.add_metabolites(metabolites)
    # scalefunc will limit the maximum and minumum fluxes
    scalefunc = lambda x: max(min(maxupperbound, x), minlowerbound)
    # read in reaction file
    reaction_file = _open_and_skip_header(baseName + ".rxn")
    reactions = []
    for line in reaction_file:
        if len(line) == 0:
            break
        the_reaction = cobra.Reaction()
        the_reaction.id = line[1]
        the_reaction.name = line[2]
        if line[3].lower() == "reversible":
            the_reaction.reversibility = 1
        elif line[3].lower() == "irreversible":
            the_reaction.reversibility = 0
        the_reaction.lower_bound = scalefunc(float(line[4]))
        the_reaction.upper_bound = scalefunc(float(line[5]))
        the_reaction.objective_coefficient = float(line[6])
        reactions.append(the_reaction)
    model.add_reactions(reactions)
    # read in S matrix
    Sfile = _open_and_skip_header(baseName + ".sto")
    S = []
    for i, line in enumerate(Sfile):
        if len(line) == 0:
            break
        the_metabolite = metabolites[i]
        for j, ns in enumerate(line):
            n = float(ns)
            if n != 0:
                model.reactions[j].add_metabolites({the_metabolite: n})
    # attempt to read in more data
    infofilepath = baseName + "_cmpd.txt"
    if maximize_info and isfile(infofilepath):
        infofile = open(infofilepath, "r")
        infofile.readline()  # skip the header
        infocsv = reader(infofile)
        for row in infocsv:
            found = _find_metabolites_by_base(row[0], model.metabolites)
            for found_metabolite in found:
                found_metabolite.formula = row[2]
                found_metabolite.parse_composition()
                found_metabolite.charge = row[4]
                found_metabolite.notes = {}
                found_metabolite.notes["KEGG_id"] = row[8]
                found_metabolite.notes["CAS"] = row[5]
                found_metaboltie.notes["review status"] = row[3]
        infofile.close()
    gpr_filepath = baseName + "_gpr.txt"
    if maximize_info and isfile(gpr_filepath):
        warn("SimPheny export files may have errors in the gpr.")
        # Using this may be risky
        gpr_file = open(gpr_filepath, "r")
        gpr_file.readline()  # skip the header
        gpr_csv = reader(gpr_file)
        for row in gpr_csv:
            the_reaction = model.reactions[model.reactions.index(row[0])]
            the_reaction.gene_reaction_rule = row[5]
            the_reaction.parse_gene_association()
        gpr_file.close()
    model.update()
    return model
