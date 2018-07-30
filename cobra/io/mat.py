# -*- coding: utf-8 -*-

from __future__ import absolute_import

import re
from collections import OrderedDict
from uuid import uuid4
from warnings import warn

from numpy import object as np_object
from numpy import array, inf, isinf
from six import string_types

from cobra.core import Metabolite, Model, Reaction
from cobra.util import create_stoichiometric_matrix
from cobra.util.solver import set_objective

try:
    import scipy.sparse as scipy_sparse
    import scipy.io as scipy_io
except ImportError:
    scipy_sparse = None
    scipy_io = None


# precompiled regular expressions
_bracket_re = re.compile(r"\[[a-z]\]$")
_underscore_re = re.compile(r"_[a-z]$")


def _get_id_compartment(id):
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
    x_no_none = [i if i is not None else "" for i in x]
    return array(x_no_none, dtype=np_object)


def load_matlab_model(infile_path, variable_name=None, inf=inf):
    """Load a cobra model stored as a .mat file

    Parameters
    ----------
    infile_path: str
        path to the file to to read
    variable_name: str, optional
        The variable name of the model in the .mat file. If this is not
        specified, then the first MATLAB variable which looks like a COBRA
        model will be used
    inf: value
        The value to use for infinite bounds. Some solvers do not handle
        infinite values so for using those, set this to a high numeric value.

    Returns
    -------
    cobra.core.Model.Model:
        The resulting cobra model

    """
    if not scipy_io:
        raise ImportError('load_matlab_model requires scipy')

    data = scipy_io.loadmat(infile_path)
    possible_names = []
    if variable_name is None:
        # skip meta variables
        meta_vars = {"__globals__", "__header__", "__version__"}
        possible_names = sorted(i for i in data if i not in meta_vars)
        if len(possible_names) == 1:
            variable_name = possible_names[0]
    if variable_name is not None:
        return from_mat_struct(data[variable_name], model_id=variable_name,
                               inf=inf)
    for possible_name in possible_names:
        try:
            return from_mat_struct(data[possible_name], model_id=possible_name,
                                   inf=inf)
        except ValueError:
            pass
    # If code here is executed, then no model was found.
    raise IOError("no COBRA model found")


def save_matlab_model(model, file_name, varname=None):
    """Save the cobra model as a .mat file.

    This .mat file can be used directly in the MATLAB version of COBRA.

    Parameters
    ----------
    model : cobra.core.Model.Model object
        The model to save
    file_name : str or file-like object
        The file to save to
    varname : string
       The name of the variable within the workspace
    """
    if not scipy_io:
        raise ImportError('load_matlab_model requires scipy')

    if varname is None:
        varname = str(model.id) \
            if model.id is not None and len(model.id) > 0 \
            else "exported_model"
    mat = create_mat_dict(model)
    scipy_io.savemat(file_name, {varname: mat},
                     appendmat=True, oned_as="column")


def create_mat_metabolite_id(model):
    for met in model.metabolites:
        if not _get_id_compartment(met.id) and met.compartment:
            yield '{}[{}]'.format(met.id,
                                  model.compartments[met.compartment].lower())
        else:
            yield met.id


def create_mat_dict(model):
    """create a dict mapping model attributes to arrays"""
    rxns = model.reactions
    mets = model.metabolites
    mat = OrderedDict()
    mat["mets"] = _cell([met_id for met_id in create_mat_metabolite_id(model)])
    mat["metNames"] = _cell(mets.list_attr("name"))
    mat["metFormulas"] = _cell([str(m.formula) for m in mets])
    try:
        mat["metCharge"] = array(mets.list_attr("charge")) * 1.
    except TypeError:
        # can't have any None entries for charge, or this will fail
        pass
    mat["genes"] = _cell(model.genes.list_attr("id"))
    # make a matrix for rxnGeneMat
    # reactions are rows, genes are columns
    rxn_gene = scipy_sparse.dok_matrix((len(model.reactions),
                                        len(model.genes)))
    if min(rxn_gene.shape) > 0:
        for i, reaction in enumerate(model.reactions):
            for gene in reaction.genes:
                rxn_gene[i, model.genes.index(gene)] = 1
        mat["rxnGeneMat"] = rxn_gene
    mat["grRules"] = _cell(rxns.list_attr("gene_reaction_rule"))
    mat["rxns"] = _cell(rxns.list_attr("id"))
    mat["rxnNames"] = _cell(rxns.list_attr("name"))
    mat["subSystems"] = _cell(rxns.list_attr("subsystem"))
    mat["csense"] = "".join((
        met._constraint_sense for met in model.metabolites))
    stoich_mat = create_stoichiometric_matrix(model)
    mat["S"] = stoich_mat if stoich_mat is not None else [[]]
    # multiply by 1 to convert to float, working around scipy bug
    # https://github.com/scipy/scipy/issues/4537
    mat["lb"] = array(rxns.list_attr("lower_bound")) * 1.
    mat["ub"] = array(rxns.list_attr("upper_bound")) * 1.
    mat["b"] = array(mets.list_attr("_bound")) * 1.
    mat["c"] = array(rxns.list_attr("objective_coefficient")) * 1.
    mat["rev"] = array(rxns.list_attr("reversibility")) * 1
    mat["description"] = str(model.id)
    return mat


def from_mat_struct(mat_struct, model_id=None, inf=inf):
    """create a model from the COBRA toolbox struct

    The struct will be a dict read in by scipy.io.loadmat

    """
    m = mat_struct
    if m.dtype.names is None:
        raise ValueError("not a valid mat struct")
    if not {"rxns", "mets", "S", "lb", "ub"} <= set(m.dtype.names):
        raise ValueError("not a valid mat struct")
    if "c" in m.dtype.names:
        c_vec = m["c"][0, 0]
    else:
        c_vec = None
        warn("objective vector 'c' not found")
    model = Model()
    if model_id is not None:
        model.id = model_id
    elif "description" in m.dtype.names:
        description = m["description"][0, 0][0]
        if not isinstance(description, string_types) and len(description) > 1:
            model.id = description[0]
            warn("Several IDs detected, only using the first.")
        else:
            model.id = description
    else:
        model.id = "imported_model"
    for i, name in enumerate(m["mets"][0, 0]):
        new_metabolite = Metabolite()
        new_metabolite.id = str(name[0][0])
        if all(var in m.dtype.names for var in
               ['metComps', 'comps', 'compNames']):
            comp_index = m["metComps"][0, 0][i][0] - 1
            new_metabolite.compartment = m['comps'][0, 0][comp_index][0][0]
            if new_metabolite.compartment not in model.compartments:
                comp_name = m['compNames'][0, 0][comp_index][0][0]
                model.compartments[new_metabolite.compartment] = comp_name
        else:
            new_metabolite.compartment = _get_id_compartment(new_metabolite.id)
            if new_metabolite.compartment not in model.compartments:
                model.compartments[
                    new_metabolite.compartment] = new_metabolite.compartment
        try:
            new_metabolite.name = str(m["metNames"][0, 0][i][0][0])
        except (IndexError, ValueError):
            pass
        try:
            new_metabolite.formula = str(m["metFormulas"][0][0][i][0][0])
        except (IndexError, ValueError):
            pass
        try:
            new_metabolite.charge = float(m["metCharge"][0, 0][i][0])
            int_charge = int(new_metabolite.charge)
            if new_metabolite.charge == int_charge:
                new_metabolite.charge = int_charge
        except (IndexError, ValueError):
            pass
        model.add_metabolites([new_metabolite])
    new_reactions = []
    coefficients = {}
    for i, name in enumerate(m["rxns"][0, 0]):
        new_reaction = Reaction()
        new_reaction.id = str(name[0][0])
        new_reaction.lower_bound = float(m["lb"][0, 0][i][0])
        new_reaction.upper_bound = float(m["ub"][0, 0][i][0])
        if isinf(new_reaction.lower_bound) and new_reaction.lower_bound < 0:
            new_reaction.lower_bound = -inf
        if isinf(new_reaction.upper_bound) and new_reaction.upper_bound > 0:
            new_reaction.upper_bound = inf
        if c_vec is not None:
            coefficients[new_reaction] = float(c_vec[i][0])
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
    set_objective(model, coefficients)
    coo = scipy_sparse.coo_matrix(m["S"][0, 0])
    for i, j, v in zip(coo.row, coo.col, coo.data):
        model.reactions[j].add_metabolites({model.metabolites[i]: v})
    return model


def _check(result):
    """ensure success of a pymatbridge operation"""
    if result["success"] is not True:
        raise RuntimeError(result["content"]["stdout"])


def model_to_pymatbridge(model, variable_name="model", matlab=None):
    """send the model to a MATLAB workspace through pymatbridge

    This model can then be manipulated through the COBRA toolbox

    Parameters
    ----------
    variable_name : str
        The variable name to which the model will be assigned in the
        MATLAB workspace

    matlab : None or pymatbridge.Matlab instance
        The MATLAB workspace to which the variable will be sent. If
        this is None, then this will be sent to the same environment
        used in IPython magics.

    """
    if scipy_sparse is None:
        raise ImportError("`model_to_pymatbridge` requires scipy!")
    if matlab is None:  # assumed to be running an IPython magic
        from IPython import get_ipython
        matlab = get_ipython().magics_manager.registry["MatlabMagics"].Matlab
    model_info = create_mat_dict(model)
    S = model_info["S"].todok()
    model_info["S"] = 0
    temp_S_name = "cobra_pymatbridge_temp_" + uuid4().hex
    _check(matlab.set_variable(variable_name, model_info))
    _check(matlab.set_variable(temp_S_name, S))
    _check(matlab.run_code("%s.S = %s;" % (variable_name, temp_S_name)))
    # all vectors need to be transposed
    for i in model_info.keys():
        if i == "S":
            continue
        _check(matlab.run_code("{0}.{1} = {0}.{1}';".format(variable_name, i)))
    _check(matlab.run_code("clear %s;" % temp_S_name))
