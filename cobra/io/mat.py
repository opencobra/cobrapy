#This section is used to load the appropriate package for writing mat files
#from Python or Jython.  Currently, only the section for Python that depends
#on scipy has been written.
from numpy import array, object as np_object
from scipy.io import loadmat, savemat
from scipy.sparse import coo_matrix


# try to use an ordered dict
try:
    from scipy.version import short_version
    scipy_version = int(short_version.split(".")[1])
    # if scipy version is earlier than 0.11, OrderedDict will not work, so use dict
    if scipy_version < 11:
        dicttype = dict
    else:
        from collections import OrderedDict as dicttype
    del short_version, scipy_version
except ImportError:
    dicttype = dict


from .. import Model, Metabolite, Reaction, Formula


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
        if m.dtype.names is None:
            continue
        if not set(["rxns", "mets", "S", "lb", "ub"]) \
                <= set(m.dtype.names):
            continue
        model = Model()
        if "description" in m:
            model.id = m["description"][0, 0][0]
        else:
            model.id = possible_name
        model.description = model.id
        for i, name in enumerate(m["mets"][0, 0]):
            new_metabolite = Metabolite()
            new_metabolite.id = str(name[0][0])
            try:
                new_metabolite.name = str(m["metNames"][0, 0][i][0][0])
                new_metabolite.formula = Formula(str(m["metFormulas"][0][0][i][0][0]))
            except:
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
                None
            try:
                new_reaction.name = str(m["rxnNames"][0, 0][i][0][0])
            except:
                pass
            new_reactions.append(new_reaction)
        model.add_reactions(new_reactions)
        coo = coo_matrix(m["S"][0, 0])
        for i, j, v in zip(coo.row, coo.col, coo.data):
            model.reactions[j].add_metabolites({model.metabolites[i]: v})
        return model
    # If code here is executed, then no model was found.
    raise Exception("no COBRA model found")


def save_matlab_model(model, file_name):
    """Save the cobra model as a .mat file.

    This .mat file can be used directly in the MATLAB version of COBRA.

    model : :class:`~cobra.core.Model.Model` object

    file_name : str or file-like object

    """
    model = model.to_array_based_model()
    rxns = model.reactions
    mets = model.metabolites
    mat = dicttype()
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
    mat["S"] = model.S
    mat["lb"] = array(rxns.list_attr("lower_bound"))
    mat["ub"] = array(rxns.list_attr("upper_bound"))
    mat["b"] = array(mets.list_attr("_bound"))
    mat["c"] = array(rxns.list_attr("objective_coefficient"))
    mat["rev"] = array(rxns.list_attr("reversibility"))
    mat["description"] = str(model.description)
    savemat(file_name, {str(model.description): mat},
             appendmat=True, oned_as="column")
