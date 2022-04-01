"""Provide functions for I/O in MATLAB (.mat) format."""

import re
from collections import OrderedDict
from typing import TYPE_CHECKING, Dict, Iterable, List, Optional
from uuid import uuid4
from warnings import warn

import numpy as np

from .. import Gene, Object
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

MET_MATLAB_TO_PROVIDERS = {
    "metHMDBID": "hmdb",
    "metInChIString": "inchi",
    "metKEGGID": "kegg.compound",
    "metPubChemID": "pubchem.compound",
    "metCHEBIID": "CHEBI",
    "metMetaNetXID": "metanetx.chemical",
    "metSEEDID": "seed.compound",
    "metBiGGID": "bigg.metabolite",
    "metBioCycID": "biocyc",
    "metEnviPathID": "envipath",
    "metLIPIDMAPSID": "lipidmaps",
    "metReactomeID": "reactome",
    "metSABIORKID": "sabiork.compound",
    "metSLMID": "SLM",
    "metSMILES": "SMILES",
    "metSBOTerm": "SBO",
}

MET_PROVIDERS_TO_MATLAB = {
    MET_MATLAB_TO_PROVIDERS[k]: k for k in MET_MATLAB_TO_PROVIDERS.keys()
}

RXN_MATLAB_TO_PROVIDERS = {
    "rxnECNumbers": "ec-code",
    "rxnReferences": "pubmed",
    "rxnKEGGID": "kegg.reaction",
    "rxnKEGGPathways": "kegg.pathway",
    "rxnMetaNetXID": "metanetx.reaction",
    "rxnSEEDID": "seed.reaction",
    "rxnBiGGID": "bigg.reaction",
    "rxnBioCycID": "biocyc",
    "rxnRheaID": "rhea",
    "rxnReactomeID": "reactome",
    "rxnSABIORKID": "sabiork.reaction",
    "rxnBRENDAID": "brenda",
    "rxnSBOTerms": "SBO",
}

RXN_PROVIDERS_TO_MATLAB = {
    RXN_MATLAB_TO_PROVIDERS[k]: k for k in RXN_MATLAB_TO_PROVIDERS.keys()
}

GENE_MATLAB_TO_PROVIDERS = {
    "geneEntrezID": "ncbigene",
    "geneRefSeqID": "refseq",
    "geneUniprotID": "uniprot",
    "geneEcoGeneID": "ecogene",
    "geneKEGGID": "kegg.gene",
    "geneHPRDID": "hprd",
    "geneASAPID": "asap",
    "geneCCDSID": "ccds",
    "geneNCBIProteinID": "ncbiprotein",
}

GENE_PROVIDERS_TO_MATLAB = {
    GENE_MATLAB_TO_PROVIDERS[k]: k for k in GENE_MATLAB_TO_PROVIDERS.keys()
}

D_GENE = "D_GENE"
D_GENE_REV = "D_GENE_REV"
D_MET = "D_MET"
D_MET_REV = "D_MET_REV"
D_REACTION = "D_REACTION"
D_REACTION_REV = "D_REACTION_REV"

D_REPLACE: dict = {
    D_GENE: GENE_PROVIDERS_TO_MATLAB,
    D_GENE_REV: GENE_MATLAB_TO_PROVIDERS,
    D_MET: MET_MATLAB_TO_PROVIDERS,
    D_MET_REV: MET_PROVIDERS_TO_MATLAB,
    D_REACTION: RXN_MATLAB_TO_PROVIDERS,
    D_REACTION_REV: RXN_PROVIDERS_TO_MATLAB,
}

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


def _cell_to_str_list(m_cell: np.ndarray, empty_value: Optional[str] = None) -> List:
    """Turn an ndarray (cell) to a list of strings.

    Parameters
    ----------
    m_cell: np.ndarray
    empty_value: str, optional
        What value to replace empty cells with. Default None.

    Returns
    -------
    List
        A list of processed strings.
    """
    return [str(x[0][0]).strip() if np.size(x[0]) else empty_value for x in m_cell]


def _cell_to_float_list(
    m_cell: np.ndarray, empty_value: Optional[float] = None
) -> List:
    """Turn an ndarray (cell) to a list of floats.

    Parameters
    ----------
    m_cell: np.ndarray
    empty_value: float, optional
        What value to replace empty cells with. Default None.

    Returns
    -------
    List
        A list of processed floats.
    """
    return [float(x[0]) if np.size(x[0]) else empty_value for x in m_cell]


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


def mat_annotations(
    target_list: List[Object], mat_struct: np.ndarray, _d_replace: str = D_MET
) -> None:
    """Process mat structure annotations in place.

    Will process mat structured annotations and add them to a list of new entities
    (metabolites, reactions, genes) in a format based on identifiers.org.

    Parameters
    ----------
    target_list: list[cobra.Object]
        A list of cobra objects, including metabolites, reactions or genes. The
        annotations will be added to these lists.
    mat_struct: np.ndarray
        A darray that includes the data imported from matlab file.
    _d_replace: str
        A string that points to the dictionary of converstions between MATLAB and
        providers. Default D_MET (for metabolite).
    """
    annotation_matlab = list(
        set(mat_struct.dtype.names).intersection(D_REPLACE[_d_replace].keys())
    )
    providers = [D_REPLACE[_d_replace][x] for x in annotation_matlab]
    annotations = dict.fromkeys(providers, None)
    for name, mat_key in zip(providers, annotation_matlab):
        annotations[name] = _cell_to_str_list(mat_struct[mat_key][0, 0])
        if mat_key == "rxnReferences":
            # noinspection PyUnresolvedReferences
            annotations[name] = [
                x.replace("PMID:", "") if x else None for x in annotations[name]
            ]
        elif mat_key == "rxnECNumbers":
            # if there are more than one ec code, turn them to a list
            # noinspection PyUnresolvedReferences,PyTypeChecker
            annotations[name] = [
                ", ".join([y.strip() for y in x.split("or")])
                if x and "or" in x
                else None
                for x in annotations[name]
            ]
        elif mat_key == "metCHEBIID":
            # if there are more than one ec code, turn them to a comma separated str
            # noinspection PyUnresolvedReferences,PyTypeChecker
            annotations[name] = [
                x.replace("CHEBI:", "") if x else None for x in annotations[name]
            ]
    for i, obj in enumerate(target_list):
        obj.annotation = {prov: annotations[prov][i]
                          for prov in providers
                          if annotations[prov][i]}


def annotations_to_mat(
    mat_struct: np.ndarray, annotation_list: List[Dict], _d_replace: str = D_MET_REV
) -> None:
    """Process mat structure annotations in place.

    Will process mat structured annotations and add them to a list of new entities
    (metabolites, reactions, genes) in a format based on identifiers.org.

    Parameters
    ----------
    mat_struct: np.ndarray
        A darray that includes the data imported from matlab file.
    annotation_list: list[Dict]
        A list of cobra annotations, in the form of a dictionary.
    _d_replace: str
        A string that points to the dictionary of converstions between MATLAB and
        providers. Default D_MET (for metabolite).
    """
    providers_used = set()
    for i in range(len(annotation_list)):
        if annotation_list[i]:
            providers_used.update(annotation_list[i].keys())
    providers_used = list(providers_used)
    annotation_matlab = [D_REPLACE[_d_replace][x] for x in providers_used]
    annotation_cells_to_be = [[None] * len(annotation_matlab) for x in annotation_list]
    for i in range(len(annotation_list)):
        if annotation_list[i]:
            for provider_key, v in annotation_list[i].items():
                provider_ind = providers_used.index(provider_key)
                annotation_cells_to_be[i][provider_ind] = v
        annotation_cells_to_be[i] = _cell(mat_struct[annotation_matlab[i]][0, 0])
        if annotation_matlab[i] == "rxnReferences":
            # noinspection PyUnresolvedReferences
            annotation_cells_to_be[i] = [
                x.replace("PMID:", "") if x else None for x in annotation_cells_to_be[i]
            ]
        if annotation_matlab[i] == "rxnECNumbers":
            # if there are more than one ec code, turn them to a comma separated str
            # noinspection PyUnresolvedReferences,PyTypeChecker
            annotation_cells_to_be[i] = [
                ", ".join([y.strip() for y in x.split("or")])
                if x and "or" in x
                else None
                for x in annotation_cells_to_be[i]
            ]
        if annotation_matlab[i] == "metCHEBIID":
            # if there are more than one ec code, turn them to a comma separated str
            # noinspection PyUnresolvedReferences,PyTypeChecker
            annotation_cells_to_be[i] = [
                x.replace("CHEBI:", "") if x else None for x in annotation_cells_to_be[i]
            ]


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
    mat_struct: np.ndarray,
    model_id: Optional[str] = None,
    inf: float = np.inf,
) -> Model:
    """Create a model from the cobratoolbox struct.

    Parameters
    ----------
    mat_struct : numpy.ndarray
        The `numpy.ndarray` that most likely contains the model, being chosen by
        load_matlab_file after loading the matlab structure via scipy_io.loadmat.
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

    if "metCharge" in m.dtype.names and "metCharges" not in m.dtype.names:
        warn(
            "This model seems to have metCharge instead of metCharges field. Will use"
            "metCharge for metabolite charges."
        )
        new_names = list(m.dtype.names)
        new_names[new_names.index("metCharge")] = "metCharges"
        m.dtype.names = new_names

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
    if "modelName" in m.dtype.names:
        model.name = m["modelName"][0, 0][0]

    met_ids = _cell_to_str_list(m["mets"][0, 0])
    if all(var in m.dtype.names for var in ["metComps", "comps", "compNames"]):
        met_comp_index = [x[0] - 1 for x in m["metComps"][0][0]]
        comps = _cell_to_str_list(m["comps"][0, 0])
        comp_names = _cell_to_str_list(m["compNames"][0][0])
        met_comps = [comps[i] for i in met_comp_index]
        met_comp_names = [comp_names[i] for i in met_comp_index]
    else:
        met_comps = [_get_id_compartment(x) for x in met_ids]
        met_comp_names = met_comps
    model.compartments.update(dict(zip(met_comps, met_comp_names)))
    met_names, met_formulas, met_charges = None, None, None
    try:
        met_names = _cell_to_str_list(m["metNames"][0, 0])
    except (IndexError, ValueError):
        # TODO: use custom cobra exception to handle exception
        pass
    try:
        met_formulas = _cell_to_str_list(m["metFormulas"][0, 0])
    except (IndexError, ValueError):
        # TODO: use custom cobra exception to handle exception
        pass
    try:
        met_charges = _cell_to_float_list(m["metCharges"][0, 0])
    except (IndexError, ValueError):
        # TODO: use custom cobra exception to handle exception
        pass
    new_metabolites = list()
    for i in range(len(met_ids)):
        new_metabolite = Metabolite(met_ids[i], compartment=met_comps[i])
        if met_names:
            new_metabolite.name = met_names[i]
        if met_charges:
            new_metabolite.charge = met_charges[i]
        if met_formulas:
            new_metabolite.formula = met_formulas[i]
        new_metabolites.append(new_metabolite)
    mat_annotations(new_metabolites, m, _d_replace=D_MET)
    model.add_metabolites(new_metabolites)

    new_genes = []
    gene_names = None
    gene_ids = _cell_to_str_list(m["genes"][0, 0])
    try:
        gene_names = _cell_to_str_list(m["geneNames"][0, 0])
    except (IndexError, ValueError):
        # TODO: use custom cobra exception to handle exception
        pass
    for i in range(len(gene_ids)):
        if gene_names:
            new_gene = Gene(gene_ids[i], name=gene_names[i])
        else:
            new_gene = Gene(gene_ids[i])
        new_genes.append(new_gene)
    mat_annotations(new_genes, m, _d_replace=D_GENE)

    new_reactions = []
    rxn_ids = _cell_to_str_list(m["rxns"][0, 0])
    rxn_lbs = _cell_to_float_list(m["lb"][0, 0])
    rxn_lbs = [-inf if np.isinf(x) and x < 0 else x for x in rxn_lbs]
    rxn_ubs = _cell_to_float_list(m["ub"][0, 0])
    rxn_ubs = [-inf if np.isinf(x) and x > 0 else x for x in rxn_ubs]
    rxn_gene_rules, rxn_names, rxn_subsystems = None, None, None
    try:
        rxn_gene_rules = _cell_to_str_list(m["grRules"][0, 0], "")
    except (IndexError, ValueError):
        # TODO: use custom cobra exception to handle exception
        pass
    try:
        rxn_names = _cell_to_str_list(m["rxnNames"][0, 0], "")
    except (IndexError, ValueError):
        # TODO: use custom cobra exception to handle exception
        pass
    try:
        rxn_subsystems = _cell_to_str_list(m["subSystems"][0, 0])
    except (IndexError, ValueError):
        # TODO: use custom cobra exception to handle exception
        pass
    for i in range(len(rxn_ids)):
        new_reaction = Reaction(
            id=rxn_ids[i], lower_bound=rxn_lbs[i], upper_bound=rxn_ubs[i]
        )
        if rxn_names:
            new_reaction.name = rxn_names[i]
        if rxn_subsystems:
            new_reaction.subsystem = rxn_subsystems[i]
        if rxn_gene_rules:
            new_reaction.gene_reaction_rule = rxn_gene_rules[i]
        new_reactions.append(new_reaction)
    mat_annotations(new_reactions, m, _d_replace=D_REACTION)

    csc = scipy_sparse.csc_matrix(m["S"][0, 0])
    for i in range(csc.shape[1]):
        stoic_dict = {
            model.metabolites[j]: csc[j, i] for j in csc.getcol(i).nonzero()[0]
        }
        new_reactions[i].add_metabolites(stoic_dict)

    model.add_reactions(new_reactions)

    if "c" in m.dtype.names:
        c_vec = _cell_to_float_list(m["c"][0, 0])
        coefficients = dict(zip(new_reactions, c_vec))
        set_objective(model, coefficients)
    else:
        warn("Objective vector `c` not found.")

    if "osenseStr" in m.dtype.names:
        model.objective_direction = str(m["osenseStr"][0, 0][0][0])
    elif "osense" in m.dtype.names:
        osense = float(m["osense"][0, 0][0][0])
        objective_direction_str = "max"
        if osense == 1:
            objective_direction_str = "min"
        model.objective_direction = objective_direction_str
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
