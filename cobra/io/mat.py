import re

from numpy import array, object as np_object
from scipy.io import loadmat, savemat
from scipy.sparse import coo_matrix

from .. import Model, Metabolite, Reaction, Formula

# try to use an ordered dict
try:
    from scipy.version import short_version
    scipy_version = int(short_version.split(".")[1])
    # if scipy version is earlier than 0.11, OrderedDict will not work, so use
    # a dict instead
    if scipy_version < 11:
        dicttype = dict
    else:
        from collections import OrderedDict as dicttype
    del short_version, scipy_version
except ImportError:
    dicttype = dict


# precompiled regular expressions
_bracket_re = re.compile("r\[[a-z]\]$")
_underscore_re = re.compile(r"_[a-z]$")


def _get_id_comparment(id):
    """extract the compartment from the id string"""
    bracket_search = _bracket_re.findall(id)
    if len(bracket_search) == 1:
        return bracket_search[0][1]
    underscore_search = _underscore_re.findall(id)
    if len(underscore_search) == 1:
        return underscore_search[0][1]
    return None


def _cell(x):
    """translate an array x into a MATLAB cell array"""
    return array(x, dtype=np_object)


def load_matlab_model(infile_path, variable_name=None):
    """Load a cobra model stored as a .mat file

    infile_path : str

    variable_name : str, optional
        The variable name of the model in the .mat file. If this is not
        specified, then the first MATLAB variable which looks like a COBRA
        model will be used

    """
    data = loadmat(infile_path)
    if variable_name is not None:
        return from_mat_struct(data[variable_name], model_id=variable_name)
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
        try:
            return from_mat_struct(data[possible_name], model_id=possible_name)
        except ValueError:
            None
    # If code here is executed, then no model was found.
    raise Exception("no COBRA model found")


def save_matlab_model(model, file_name):
    """Save the cobra model as a .mat file.

    This .mat file can be used directly in the MATLAB version of COBRA.

    model : :class:`~cobra.core.Model.Model` object

    file_name : str or file-like object

    """
    mat = create_mat_dict(model)
    savemat(file_name, {str(model.description): mat},
            appendmat=True, oned_as="column")


def create_mat_dict(model):
    """create a dict mapping model attributes to arrays"""
    model = model.to_array_based_model()
    rxns = model.reactions
    mets = model.metabolites
    mat = dicttype()
    mat["mets"] = _cell(mets.list_attr("id"))
    mat["metNames"] = _cell(mets.list_attr("name"))
    mat["metFormulas"] = _cell([str(m.formula) for m in mets])
    mat["genes"] = _cell(model.genes.list_attr("id"))
    mat["grRules"] = _cell(rxns.list_attr("gene_reaction_rule"))
    mat["rxns"] = _cell(rxns.list_attr("id"))
    mat["rxnNames"] = _cell(rxns.list_attr("name"))
    mat["subSystems"] = _cell(rxns.list_attr("subsystem"))
    mat["csense"] = "".join(model._constraint_sense)
    mat["S"] = model.S
    mat["lb"] = array(rxns.list_attr("lower_bound"))
    mat["ub"] = array(rxns.list_attr("upper_bound"))
    mat["b"] = array(mets.list_attr("_bound"))
    mat["c"] = array(rxns.list_attr("objective_coefficient"))
    mat["rev"] = array(rxns.list_attr("reversibility"))
    mat["description"] = str(model.description)
    return mat


def from_mat_struct(mat_struct, model_id=None):
    """create a model from the COBRA toolbox struct

    The struct will be a dict read in by scipy.io.loadmat

    """
    m = mat_struct
    if m.dtype.names is None:
        raise ValueError("not a valid mat struct")
    if not set(["rxns", "mets", "S", "lb", "ub"]) <= set(m.dtype.names):
        raise ValueError("not a valid mat struct")
    model = Model()
    if "description" in m:
        model.id = m["description"][0, 0][0]
    elif model_id is not None:
        model.id = model_id
    else:
        model.id = "imported_model"
    model.description = model.id
    for i, name in enumerate(m["mets"][0, 0]):
        new_metabolite = Metabolite()
        new_metabolite.id = str(name[0][0])
        new_metabolite.compartment = _get_id_comparment(new_metabolite.id)
        try:
            new_metabolite.name = str(m["metNames"][0, 0][i][0][0])
        except (IndexError, ValueError):
            pass
        try:
            new_metabolite.formula = \
                Formula(str(m["metFormulas"][0][0][i][0][0]))
        except (IndexError, ValueError):
            pass
        model.add_metabolites([new_metabolite])
    new_reactions = []
    for i, name in enumerate(m["rxns"][0, 0]):
        new_reaction = Reaction()
        new_reaction.id = str(name[0][0])
        new_reaction.lower_bound = float(m["lb"][0, 0][i][0])
        new_reaction.upper_bound = float(m["ub"][0, 0][i][0])
        new_reaction.objective_coefficient = float(m["c"][0, 0][i][0])
        try:
            new_reaction.gene_reaction_rule = str(m['grRules'][0, 0][i][0][0])
        except (IndexError, ValueError):
            pass
        try:
            new_reaction.name = str(m["rxnNames"][0, 0][i][0][0])
        except (IndexError, ValueError):
            pass
        try:
            new_reaction.subsystem = str(m['subSystems'][0, 0][i][0][0])
        except (IndexError, ValueError):
            pass
        new_reactions.append(new_reaction)
    model.add_reactions(new_reactions)
    coo = coo_matrix(m["S"][0, 0])
    for i, j, v in zip(coo.row, coo.col, coo.data):
        model.reactions[j].add_metabolites({model.metabolites[i]: v})
    return model
