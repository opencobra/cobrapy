"""Provide functions for I/O in MATLAB (.mat) format."""

import re
from collections import OrderedDict
from typing import TYPE_CHECKING, Dict, Iterable, Optional
from uuid import uuid4
from warnings import warn

import numpy as np

from ..core import Metabolite, Model, Reaction
from ..util import create_stoichiometric_matrix
from ..util.solver import set_objective


try:
    import scipy.io as scipy_io
    import scipy.sparse as scipy_sparse
except ImportError:
    scipy_sparse = None
    scipy_io = None


if TYPE_CHECKING:
    import pymatbridge


# precompiled regular expressions (kept globally for caching)
_bracket_re = re.compile(r"\[(?P<compartment>[a-z]+)\]$")
_underscore_re = re.compile(r"_(?P<compartment>[a-z]+)$")


def _get_id_compartment(id: str) -> str:
    """Extract the compartment from the `id` string.

    Parameters
    ----------
    id : str
        The ID string to extract component from.

    Returns
    -------
    str
        The extracted component string.

    """
    bracket_search = _bracket_re.search(id)
    if bracket_search:
        return bracket_search.group("compartment")

    underscore_search = _underscore_re.search(id)
    if underscore_search:
        return underscore_search.group("compartment")


def _cell(x: Iterable[str]) -> np.ndarray:
    """Translate an iterable `x` into a MATLAB cell array.

    Parameters
    ----------
    x : iterable of str
        The data iterable to convert to cell array.

    Returns
    -------
    numpy.ndarray
       The converted cell array compatible with MATLAB.

    """
    x_no_none = [i if i is not None else "" for i in x]
    return np.array(x_no_none, dtype=object)


def load_matlab_model(
    infile_path: str, variable_name: Optional[str] = None, inf: float = np.inf
) -> Model:
    """Load a cobra model stored as a .mat file.

    Parameters
    ----------
    infile_path : str
        File path or descriptor of the .mat file describing the cobra model.
    variable_name : str, optional
        The variable name of the model in the .mat file. If None, then the
        first MATLAB variable which looks like a COBRA model will be used
        (default None).
    inf: float, optional
        The value to use for infinite bounds. Some solvers do not handle
        infinite values so for using those, set this to a high numeric value
        (default `numpy.inf`).

    Returns
    -------
    cobra.Model
        The cobra model as represented in the .mat file.

    Raises
    ------
    ImportError
        If scipy is not found in the Python environment.
    IOError
        If no COBRA model is found in the .mat file.

    """
    if not scipy_io:
        raise ImportError("load_matlab_model() requires scipy.")

    data = scipy_io.loadmat(infile_path)
    possible_names = []
    if variable_name is None:
        # skip meta variables
        meta_vars = {"__globals__", "__header__", "__version__"}
        possible_names = sorted(i for i in data if i not in meta_vars)
        if len(possible_names) == 1:
            variable_name = possible_names[0]
    elif variable_name is not None:
        return from_mat_struct(data[variable_name], model_id=variable_name, inf=inf)

    for possible_name in possible_names:
        try:
            return from_mat_struct(data[possible_name], model_id=possible_name, inf=inf)
        except ValueError:
            # TODO: use custom cobra exception to handle exception
            pass
    # If code here is executed, then no model was found.
    raise IOError(f"No COBRA model found at {infile_path}.")


def save_matlab_model(
    model: Model, file_name: str, varname: Optional[str] = None
) -> None:
    """Save the cobra model as a .mat file.

    This .mat file can be used directly in cobratoolbox.

    Parameters
    ----------
    model : cobra.Model
        The cobra model to represent.
    file_name : str or file-like
        File path or descriptor that the MATLAB representation should be
        written to.
    varname : str, optional
       The name of the variable within the MATLAB workspace. Model ID is
       used if available, else 'exported_model' is used (default None).

    """
    if not scipy_io:
        raise ImportError("save_matlab_model() requires scipy.")

    if varname is None:
        varname = (
            str(model.id)
            if model.id is not None and len(model.id) > 0
            else "exported_model"
        )
    mat = create_mat_dict(model)
    scipy_io.savemat(file_name, {varname: mat}, appendmat=True, oned_as="column")


def create_mat_metabolite_id(model: Model) -> str:
    """Obtain all metabolite IDs from a MATLAB model.

    Parameters
    ----------
    model : cobra.Model
        The model to obtain metabolite IDs from.

    Yields
    ------
    str
        The metabolite ID along with compartment (if found).

    """
    for met in model.metabolites:
        if not _get_id_compartment(met.id) and met.compartment:
            yield "{}[{}]".format(met.id, model.compartments[met.compartment].lower())
        else:
            yield met.id


def create_mat_dict(model: Model) -> OrderedDict:
    """Create a dictionary mapping model attributes to arrays.

    Parameters
    ----------
    model : cobra.Model
        The model to create dictionary for.

    Returns
    -------
    OrderedDict
        The ordered dictionary having model attributes as keys and their
        respective values represented as arrays, as the values.

    """
    rxns = model.reactions
    mets = model.metabolites
    mat = OrderedDict()
    mat["mets"] = _cell([met_id for met_id in create_mat_metabolite_id(model)])
    mat["metNames"] = _cell(mets.list_attr("name"))
    mat["metFormulas"] = _cell([str(m.formula) for m in mets])
    try:
        mat["metCharge"] = np.array(mets.list_attr("charge")) * 1.0
    except (TypeError, AttributeError):
        # can't have any None entries for charge, or this will fail
        # TODO: use custom cobra exception to handle exception
        pass
    mat["genes"] = _cell(model.genes.list_attr("id"))
    # make a matrix for rxnGeneMat
    # reactions are rows, genes are columns
    rxn_gene = scipy_sparse.dok_matrix((len(model.reactions), len(model.genes)))
    if min(rxn_gene.shape) > 0:
        for i, reaction in enumerate(model.reactions):
            for gene in reaction.genes:
                rxn_gene[i, model.genes.index(gene)] = 1
        mat["rxnGeneMat"] = rxn_gene
    mat["grRules"] = _cell(rxns.list_attr("gene_reaction_rule"))
    mat["rxns"] = _cell(rxns.list_attr("id"))
    mat["rxnNames"] = _cell(rxns.list_attr("name"))
    mat["subSystems"] = _cell(rxns.list_attr("subsystem"))
    stoich_mat = create_stoichiometric_matrix(model)
    mat["S"] = stoich_mat if stoich_mat is not None else [[]]
    # multiply by 1 to convert to float, working around scipy bug
    # https://github.com/scipy/scipy/issues/4537
    mat["lb"] = np.array(rxns.list_attr("lower_bound")) * 1.0
    mat["ub"] = np.array(rxns.list_attr("upper_bound")) * 1.0
    mat["b"] = np.array(mets.list_attr("_bound")) * 1.0
    mat["c"] = np.array(rxns.list_attr("objective_coefficient")) * 1.0
    mat["rev"] = np.array(rxns.list_attr("reversibility")) * 1
    mat["description"] = str(model.id)
    return mat


def from_mat_struct(
    mat_struct: Dict[str, np.ndarray],
    model_id: Optional[str] = None,
    inf: float = np.inf,
) -> Model:
    """Create a model from the cobratoolbox struct.

    Parameters
    ----------
    mat_struct : dict
        The dictionary loaded via `scipy.io.loadmat`, having str as keys
        and `numpy.ndarray` as values.
    model_id : str, optional
        The ID of the model generated. If None, will try to look for ID in
        model's description. If multiple IDs are found, the first one is
        used. If no IDs are found, will use 'imported_model' (default None).
    inf : float, optional
        The value to use for infinite bounds. Some solvers do not handle
        infinite values so for using those, set this to a high numeric value
        (default `numpy.inf`).

    Returns
    -------
    cobra.Model
        The model as represented in .mat file.

    """
    m = mat_struct

    if m.dtype.names is None or not {"rxns", "mets", "S", "lb", "ub"} <= set(
        m.dtype.names
    ):
        raise ValueError("Invalid MATLAB struct.")

    if "c" in m.dtype.names:
        c_vec = m["c"][0, 0]
    else:
        c_vec = None
        warn("Objective vector `c` not found.")

    model = Model()
    if model_id is not None:
        model.id = model_id
    elif "description" in m.dtype.names:
        description = m["description"][0, 0][0]
        if not isinstance(description, str) and len(description) > 1:
            model.id = description[0]
            warn("Several IDs detected, only using the first.")
        else:
            model.id = description
    else:
        model.id = "imported_model"

    new_metabolites = list()
    for i, name in enumerate(m["mets"][0, 0]):
        new_metabolite = Metabolite()
        new_metabolite.id = str(name[0][0]).strip()
        if all(var in m.dtype.names for var in ["metComps", "comps", "compNames"]):
            comp_index = m["metComps"][0, 0][i][0] - 1
            new_metabolite.compartment = m["comps"][0, 0][comp_index][0][0]
            if new_metabolite.compartment not in model.compartments:
                comp_name = m["compNames"][0, 0][comp_index][0][0]
                model.compartments[new_metabolite.compartment] = comp_name
        else:
            new_metabolite.compartment = _get_id_compartment(new_metabolite.id)
            if new_metabolite.compartment not in model.compartments:
                model.compartments[
                    new_metabolite.compartment
                ] = new_metabolite.compartment
        try:
            new_metabolite.name = str(m["metNames"][0, 0][i][0][0])
        except (IndexError, ValueError):
            # TODO: use custom cobra exception to handle exception
            pass
        try:
            new_metabolite.formula = str(m["metFormulas"][0][0][i][0][0])
        except (IndexError, ValueError):
            # TODO: use custom cobra exception to handle exception
            pass
        try:
            new_metabolite.charge = float(m["metCharge"][0, 0][i][0])
            int_charge = int(new_metabolite.charge)
            if new_metabolite.charge == int_charge:
                new_metabolite.charge = int_charge
        except (IndexError, ValueError):
            # TODO: use custom cobra exception to handle exception
            pass
        new_metabolites.append(new_metabolite)
    model.add_metabolites(new_metabolites)

    new_reactions = []
    coefficients = {}
    for i, name in enumerate(m["rxns"][0, 0]):
        new_reaction = Reaction()
        new_reaction.id = str(name[0][0]).strip()
        new_reaction.bounds = float(m["lb"][0, 0][i][0]), float(m["ub"][0, 0][i][0])
        if np.isinf(new_reaction.lower_bound) and new_reaction.lower_bound < 0:
            new_reaction.lower_bound = -inf
        if np.isinf(new_reaction.upper_bound) and new_reaction.upper_bound > 0:
            new_reaction.upper_bound = inf
        if c_vec is not None:
            coefficients[new_reaction] = float(c_vec[i][0])
        try:
            new_reaction.gene_reaction_rule = str(m["grRules"][0, 0][i][0][0]).strip()
        except (IndexError, ValueError):
            # TODO: use custom cobra exception to handle exception
            pass
        try:
            new_reaction.name = str(m["rxnNames"][0, 0][i][0][0]).strip()
        except (IndexError, ValueError):
            # TODO: use custom cobra exception to handle exception
            pass
        try:
            new_reaction.subsystem = str(m["subSystems"][0, 0][i][0][0]).strip()
        except (IndexError, ValueError):
            # TODO: use custom cobra exception to handle exception
            pass
        new_reactions.append(new_reaction)
    model.add_reactions(new_reactions)
    set_objective(model, coefficients)

    csc = scipy_sparse.csc_matrix(m["S"][0, 0])
    for i in range(csc.shape[1]):
        stoic_dict = {model.metabolites[j]: csc[j, i] for j in
                      csc.getcol(i).nonzero()[0]}
        model.reactions[i].add_metabolites(stoic_dict)
    return model


def _check(result: Dict[str, str]) -> None:
    """Ensure success of a `pymatbridge` operation.

    Parameters
    ----------
    result : dict
        The dictionary obtained from `pymatbridge` with keys as 'message',
        'result', and 'success'.

    Raises
    ------
    RuntimeError
        If setting variable via `pymatbridge` fails.

    """
    if result["success"] is not True:
        # TODO: verify if key 'content' is valid as docs don't state about it
        raise RuntimeError(result["content"]["stdout"])


def model_to_pymatbridge(
    model: Model,
    variable_name: str = "model",
    matlab: Optional["pymatbridge.Matlab"] = None,
) -> None:
    """Send the model to a MATLAB workspace through `pymatbridge`.

    This model can then be manipulated through the cobratoolbox.

    Parameters
    ----------
    model : cobra.Model
        The model to send to MATLAB workspace.
    variable_name : str, optional
        The variable name to which the model will be assigned in the
        MATLAB workspace (default 'model').
    matlab : pymatbridge.Matlab, optional
        The MATLAB workspace to which the variable will be sent. If None,
        the variable will be sent to the same environment used in IPython
        magics.

    """
    if scipy_sparse is None:
        raise ImportError("model_to_pymatbridge() requires scipy.")

    if matlab is None:  # assumed to be running an IPython magic
        from IPython import get_ipython

        matlab = get_ipython().magics_manager.registry["MatlabMagics"].Matlab

    model_info = create_mat_dict(model)
    S = scipy_sparse.dok_matrix(model_info["S"])
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
