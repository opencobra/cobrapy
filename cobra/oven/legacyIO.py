from os.path import isfile
import csv
import re
from warnings import warn
try:
    from cPickle import load, dump
except ImportError:
    from pickle import load, dump
# Because OrderedDict is not implemented on all systems
try:
    from collections import OrderedDict as _dicttype
except ImportError:
    _dicttype = dict
from numpy import array as _array, object as _object
from scipy.io import loadmat as _loadmat, savemat as _savemat
from scipy.sparse import coo_matrix as _coo_matrix
import cobra
from cobra.io.sbml import create_cobra_model_from_sbml_file as _sbml_import


# function to help translate an array to a MATLAB cell array
_cell = lambda x: _array(x, dtype=_object)


def load_pickle(pickle_file):
    """Read in a pickle file.

    Parameters
    ----------
    pickle_file : str
        The path to the pickle to load.

    Returns
    -------
    contents : the contents of the pickle

    """
    # if the user does not add the .pickle extension
    if not isfile(pickle_file) and not pickle_file.endswith(".pickle") \
            and isfile("pickle_file" + ".pickle"):
        pickle_file += ".pickle"
    infile = open(pickle_file, "rb")
    contents = load(infile)
    infile.close()
    return contents


def export_flux_distribution(model, filepath):
    """Export flux distribution to import into Simpheny.

    Parameters
    ----------
    model : cobra.Model
    filepath: str

    """
    from simphenyMapping import mapping
    outfile = open(filepath, "w")
    outcsv = csv.writer(outfile, delimiter="\t", lineterminator="\n")
    outcsv.writerow(["Reaction Number", "Flux Value",
                     "Lower Bound", "Upper Bound"])
    for reaction_name, reaction_flux in model.solution.x_dict.iteritems():
        reaction = model.reactions.get_by_id(reaction_name)
        try:
            outcsv.writerow([mapping[reaction_name], reaction_flux,
                reaction.lower_bound, reaction.upper_bound])
        except KeyError, e:
            print "Simpheny id number not found for", e
    outfile.close()


def load_matlab_model(infile_path, variable_name=None):
    """Load a cobra model stored as a .mat file
    NOTE: INCOMPLETE, does not load GPR's

    Parameters
    ----------
    infile_path : str
    variable_name : str, optional
        The variable name of the model in the .mat file. If this is not
        specified, then the first MATLAB variable which looks like a COBRA
        model will be used

    """
    data = _loadmat(infile_path)
    if variable_name is not None:
        possible_names = [variable_name]
    else:
        # will try all of the variables in the dict
        possible_names = {}
        for key in data.keys():
            possible_names[key] = None
        # skip meta variables
        to_remove = ["__globals__", "__header__", "__version__"]
        to_pop = []
        for name in possible_names:
            if name in to_remove:
                to_pop.append(name)
        for i in to_pop:
            possible_names.pop(i)
        possible_names = possible_names.keys()
    for possible_name in possible_names:
        m = data[possible_name]  # TODO: generalize
        if not set(["rxns", "mets", "S", "lb", "ub"]) \
                <= set(m.dtype.names):
            continue
        model = cobra.Model()
        model.id = m["description"][0, 0][0]
        model.description = model.id
        for i, name in enumerate(m["mets"][0, 0]):
            new_metabolite = cobra.Metabolite()
            new_metabolite.id = name[0][0]
            new_metabolite.name = m["metNames"][0, 0][i][0][0]
            new_metabolite.formula = m["metFormulas"][0][0][i][0][0]
            model.metabolites.append(new_metabolite)
        for i, name in enumerate(m["rxns"][0, 0]):
            new_reaction = cobra.Reaction()
            new_reaction.id = name[0][0]
            new_reaction.name = m["rxnNames"][0, 0][i][0][0]
            new_reaction.lower_bound = m["lb"][0, 0][i][0]
            new_reaction.upper_bound = m["ub"][0, 0][i][0]
            new_reaction.objective_coefficient = m["c"][0, 0][i][0]
            model.reactions.append(new_reaction)
        coo = _coo_matrix(m["S"][0, 0])
        for i, j, v in zip(coo.row, coo.col, coo.data):
            model.reactions[j].add_metabolites({model.metabolites[i]: v})
        # TODO finish adding GPR's
        # model.update()
        return model
    # If code here is executed, then no model was found.
    raise Exception("no COBRA model found")


def save_matlab_model(model, outfile_path):
    """Save the cobra model as a .mat file.

    This .mat file can be used directly in the MATLAB version of COBRA.
    NOTE: using this function requires a patched version of cobra.io.savemat

    Parameters
    ----------
    model : cobra.Model
    outfile_path : str

    """
    model.update()
    rxns = model.reactions
    mets = model.metabolites
    mat = _dicttype()
    csense = ""
    for m in mets:
        csense += m._constraint_sense
    mat["mets"] = _cell(mets.list_attr("id"))
    mat["metNames"] = _cell(mets.list_attr("name"))
    mat["metFormulas"] = _cell([str(m.formula) for m in mets])
    mat["genes"] = _cell(model.genes.list_attr("id"))
    mat["grRules"] = _cell(rxns.list_attr("gene_reaction_rule"))
    mat["rxns"] = _cell(rxns.list_attr("id"))
    mat["rxnNames"] = _cell(rxns.list_attr("name"))
    mat["subSystems"] = _cell(rxns.list_attr("subsystem"))
    mat["csense"] = csense
    mat["S"] = model._S
    mat["lb"] = _array(rxns.list_attr("lower_bound"))
    mat["ub"] = _array(rxns.list_attr("upper_bound"))
    mat["b"] = _array(mets.list_attr("_bound"))
    mat["c"] = _array(rxns.list_attr("objective_coefficient"))
    mat["rev"] = _array(rxns.list_attr("reversibility"))
    mat["description"] = str(model.description)
    _savemat(outfile_path, {str(model.description): mat},
             appendmat=True, oned_as="column")


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
    the_id = the_id.replace('__', '-')
    return the_id


def read_legacy_sbml(filename):
    """read in an sbml file and fix the sbml id's"""
    model = _sbml_import(filename)
    for metabolite in model.metabolites:
        metabolite.id = _fix_legacy_id(metabolite.id)
    model.metabolites._generate_index()
    for reaction in model.reactions:
        reaction.id = _fix_legacy_id(reaction.id)
        reaction.reconstruct_reaction()
    model.reactions._generate_index()
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
    return csv.reader(file, delimiter="\t")


def _find_metabolites_by_base(base, metabolites):
    """search for all metabolites in the list which match the base name.
    For example, "h2o" will identify both "h2o(c)" and "h2o(e)" """
    search = re.compile(base + "\([a-z]\)")
    found = []
    for the_metabolite in metabolites:
        if search.match(the_metabolite.id) is not None:
            found.append(the_metabolite)
    return found


def read_simpheny(baseName, min_lower_bound=-1000, max_upper_bound=1000,
        maximize_info=True):
    r"""Imports files exported from a SimPheny simulation as a cobra model.

    .. warning:: Use with caution. This is a legacy import function, and
        errors have been observed in the converted gene-reaction rules.

    Parameters
    ----------
    baseName : str
        The filepath to the exported SimPheny files without any of the
        extensions. On Windows, it helps if baseName is a raw string
        (i.e. r"Path\\to\\files")
    min_lower_bound, max_upper_bound : float or int, optional
        The bounds on the lower and upper bounds of fluxes in model.
    maximize_info : bool, optional
        An optional boolean keyword argument. If True, then an attempt
        will be made to parse the gpr and metabolite info files, and the
        function will take a little bit longer.

    Returns
    -------
    model : cobra.Model
        the imported simpheny model

    """

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
            compartment_search = re.findall("\([a-z]\)$", metabolite.id)
            if compartment_search != []:
                metabolite.compartment = compartment_search[0][1]
                model.compartments[metabolite.compartment] = line[3]
        metabolites.append(metabolite)
    model.add_metabolites(metabolites)
    # scalefunc will limit the maximum and minumum fluxes
    scalefunc = lambda x: max(min(max_upper_bound, x), min_lower_bound)
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
        infocsv = csv.reader(infofile)
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
        gpr_csv = csv.reader(gpr_file)
        for row in gpr_csv:
            the_reaction = model.reactions[model.reactions.index(row[0])]
            the_reaction.gene_reaction_rule = row[5]
            the_reaction.parse_gene_association()
        gpr_file.close()
    # model.update()
    return model

if __name__ == "__main__":
    model = load_matlab_model(r"C:\Users\aebrahim\Documents\MATLAB\iJO1366")
    model.update()
    # from cobra.test import ecoli_pickle
    # model = load_pickle(ecoli_pickle)
    # save_matlab_model(model, r"C:\Users\aebrahim\Documents\Matlab\test.mat")
