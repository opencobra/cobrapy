"""Provide functions for I/O in MATLAB (.mat) format."""

import logging
import re
from collections import OrderedDict
from pathlib import Path
from typing import IO, Dict, Iterable, List, Optional, Pattern, Union

import numpy as np

from ..core import Gene, Group, Metabolite, Model, Object, Reaction
from ..util import create_stoichiometric_matrix
from ..util.solver import set_objective


try:
    import scipy.io as scipy_io
    import scipy.sparse as scipy_sparse
except ImportError:
    scipy_sparse = None
    scipy_io = None

logger = logging.getLogger(__name__)

# The following dictionaries are based on
# https://github.com/opencobra/cobratoolbox/blob/docs/source/notes/COBRAModelFields.md
# at commit 83d26938a9babff79289d40e20f5f50dd5b710fa
# Which is the most updated master as of April 4th, 2022

MET_MATLAB_TO_PROVIDERS = {
    "metHMDBID": "hmdb",
    "metInChIString": "inchi",
    "metKEGGID": "kegg.compound",
    "metKEGGGlycanID": "kegg.glycan",
    "metKEGGDrugID": "kegg.drug",
    "metUniPathway": "unipathway.compound",
    "metPubChemID": "pubchem.compound",
    "metPubChemSubstance": "pubchem.substance",
    "metCHEBIID": "chebi",
    "metMetaNetXID": "metanetx.chemical",
    "metSEEDID": "seed.compound",
    "metBiGGID": "bigg.metabolite",
    "metBioCycID": "biocyc",
    "metEnviPathID": "envipath",
    "metLIPIDMAPSID": "lipidmaps",
    "metReactomeID": "reactome",
    "metSABIORKID": "sabiork.compound",
    "metSLMID": "slm",
    "metSMILES": "SMILES",
    "metSBOTerms": "sbo",
    "metCasNumber": "cas",
}

MET_PROVIDERS_TO_MATLAB = {
    MET_MATLAB_TO_PROVIDERS[k]: k for k in MET_MATLAB_TO_PROVIDERS.keys()
}

MET_NOTES_TO_MATLAB = {
    "metNotes": "metNotes",
}

MET_MATLAB_TO_NOTES = {MET_NOTES_TO_MATLAB[k]: k for k in MET_NOTES_TO_MATLAB.keys()}

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
    "rxnSBOTerms": "sbo",
}

RXN_PROVIDERS_TO_MATLAB = {
    RXN_MATLAB_TO_PROVIDERS[k]: k for k in RXN_MATLAB_TO_PROVIDERS.keys()
}

CONFIDENCE_STR = "Confidence Level"
RXN_MATLAB_TO_NOTES = {
    "rxnReferences": "References",
    "rxnNotes": "NOTES",
    "rxnConfidenceScores": CONFIDENCE_STR,
}

RXN_NOTES_TO_MATLAB = {
    CONFIDENCE_STR: "rxnConfidenceScores",
    "NOTES": "rxnNotes",
    "References": "rxnNotes",
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

DICT_GENE = "DICT_GENE"
DICT_GENE_REV = "DICT_GENE_REV"
DICT_MET = "DICT_MET"
DICT_MET_REV = "DICT_MET_REV"
DICT_MET_NOTES = "DICT_MET_NOTES"
DICT_MET_NOTES_REV = "DICT_MET_NOTES_REV"
DICT_REACTION = "DICT_REACTION"
DICT_REACTION_REV = "DICT_REACTION_REV"
DICT_REACTION_NOTES = "DICT_REACTION_NOTES"
DICT_REACTION_NOTES_REV = "DICT_REACTION_NOTES_REV"

DICT_REPLACE: dict = {
    DICT_GENE: GENE_MATLAB_TO_PROVIDERS,
    DICT_GENE_REV: GENE_PROVIDERS_TO_MATLAB,
    DICT_MET: MET_MATLAB_TO_PROVIDERS,
    DICT_MET_REV: MET_PROVIDERS_TO_MATLAB,
    DICT_MET_NOTES: MET_MATLAB_TO_NOTES,
    DICT_MET_NOTES_REV: MET_NOTES_TO_MATLAB,
    DICT_REACTION: RXN_MATLAB_TO_PROVIDERS,
    DICT_REACTION_REV: RXN_PROVIDERS_TO_MATLAB,
    DICT_REACTION_NOTES: RXN_MATLAB_TO_NOTES,
    DICT_REACTION_NOTES_REV: RXN_NOTES_TO_MATLAB,
}

# precompiled regular expressions (kept globally for caching)
_bracket_re = re.compile(r"[\[(](?P<compartment>[a-zA-Z]+)[])]$")
# Some older models have _boundary for (some) exchange reactions
_underscore_re = re.compile(r"_(?P<compartment>[a-zA-Z]+)(_boundary)?$")
_pubmed_re = re.compile("PMID: ?(\\d+),?")
_punctuation_re = re.compile(r"^[;,.'\"]+$")
_double_punctuation_re = re.compile(r"[;,.'\"]{2,}")
_ec_re = re.compile(r"([\d\-]+.[\d\-]+.[\d\-]+.[\d-]+)")
_chebi_re = re.compile(r"\D*(\d+),?")
_sbo_re = re.compile(r"\D*(\d+),?")


def _get_id_compartment(_id: str) -> str:
    """Extract the compartment from the `id` string.

    Parameters
    ----------
    _id : str
        The ID string to extract component from.

    Returns
    -------
    str
        The extracted component string.

    """
    bracket_search = _bracket_re.search(_id)
    if bracket_search:
        return bracket_search.group("compartment")

    underscore_search = _underscore_re.search(_id)
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


def _cell_to_str_list(
    m_cell: np.ndarray,
    empty_value: Optional[str] = None,
    pattern_split: Optional[Pattern] = None,
    str_prefix: str = "",
) -> List:
    """Turn an ndarray (cell) to a list of strings.

    Parameters
    ----------
    m_cell: np.ndarray
    empty_value: str, optional
        What value to replace empty cells with. Default None.
    pattern_split: Pattern, optional
        Regular expression to use to split the expression. Used for annotations.
        Default None.
    str_prefix: str, optional
        A prefix that will be added to each value in the list if present. Default "".

    Returns
    -------
    List
        A list of processed strings.
    """
    if str_prefix and pattern_split:
        return [
            [
                str_prefix + str_found if str_prefix not in str_found else str_found
                for str_found in pattern_split.findall(str(each_cell[0][0]))
            ]
            if np.size(each_cell[0])
            else empty_value
            for each_cell in m_cell
        ]
    elif pattern_split:
        return [
            pattern_split.findall(str(each_cell[0][0]))
            if np.size(each_cell[0])
            else empty_value
            for each_cell in m_cell
        ]
    else:
        return [
            str(each_cell[0][0]).strip() if np.size(each_cell[0]) else empty_value
            for each_cell in m_cell
        ]


def _cell_to_float_list(
    m_cell: np.ndarray,
    empty_value: Optional[float] = None,
    inf_value: Optional[float] = None,
) -> List:
    """Turn an ndarray (cell) to a list of floats.

    Parameters
    ----------
    m_cell: np.ndarray
    empty_value: float, optional
        What value to replace empty cells with. Default None.
    inf_value: float, optional
        Replace infinite values with defined inf. Default None (don't replace), will be
        replaced by given value  if given.

    Returns
    -------
    List
        A list of processed floats.
    """

    def fix_inf(val: str, _inf_value: float):
        """Fix inf value, used for rxn.lb and rxn.ub.

        Parameters
        ----------
        val: str
            the str to process into float
        _inf_value: float
            A value to replace infite values with.
        """
        val = float(val)
        if np.isinf(val) and val < 0:
            return -_inf_value
        elif np.isinf(val) and val > 0:
            return _inf_value
        else:
            return val

    if inf_value:
        return [
            fix_inf(x[0], _inf_value=inf_value) if np.size(x[0]) else empty_value
            for x in m_cell
        ]
    else:
        return [float(x[0]) if np.size(x[0]) else empty_value for x in m_cell]


def load_matlab_model(
    infile_path: Union[str, Path, IO],
    variable_name: Optional[str] = None,
    inf: float = np.inf,
) -> Model:
    """Load a cobra model stored as a .mat file.

    Parameters
    ----------
    infile_path : str or Path or filehandle
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

    if isinstance(infile_path, str):
        data = scipy_io.loadmat(infile_path)
    elif isinstance(infile_path, Path):
        data = scipy_io.loadmat(infile_path.open("rb"))  # noqa W9018
    else:
        data = scipy_io.loadmat(infile_path)  # noqa W9018
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
        except ValueError as e:
            print(f"Some problem with the model, causing error {e}")
            # TODO: use custom cobra exception to handle exception
    # If code here is executed, then no model was found.
    raise IOError(f"No COBRA model found at {infile_path}.")


def save_matlab_model(
    model: Model, file_name: Union[str, Path, IO], varname: Optional[str] = None
) -> None:
    """Save the cobra model as a .mat file.

    This .mat file can be used directly in cobratoolbox.

    Parameters
    ----------
    model : cobra.Model
        The cobra model to represent.
    file_name : str or file-like or Path
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


def mat_parse_annotations(
    target_list: List[Object], mat_struct: np.ndarray, d_replace: str = DICT_MET
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
    d_replace: str
        A string that points to the dictionary of converstions between MATLAB and
        providers. Default DICT_MET (for metabolite).
    """
    struct_names = [x.casefold() for x in mat_struct.dtype.names]
    matlab_field_dict = {
        x.casefold(): DICT_REPLACE[d_replace][x] for x in DICT_REPLACE[d_replace].keys()
    }
    caseunfold = {x.casefold(): x for x in mat_struct.dtype.names}
    annotation_matlab = list(set(struct_names).intersection(matlab_field_dict.keys()))
    providers = [matlab_field_dict[x] for x in annotation_matlab]
    annotations = dict.fromkeys(providers, None)
    for name, mat_key in zip(providers, annotation_matlab):
        if mat_key == "rxnReferences".casefold():
            # This only picks up PMID: style references. Sometimes there are other
            # things like PMC or OMIM, but those are ignored for now,
            annotations[name] = _cell_to_str_list(
                mat_struct[caseunfold[mat_key]][0, 0], None, _pubmed_re
            )
        elif mat_key == "rxnECNumbers".casefold():
            # turn EC codes to a list
            annotations[name] = _cell_to_str_list(
                mat_struct[caseunfold[mat_key]][0, 0], None, _ec_re
            )
        elif mat_key == "metCHEBIID".casefold():
            annotations[name] = _cell_to_str_list(
                mat_struct[caseunfold[mat_key]][0, 0], None, _chebi_re, "CHEBI:"
            )
        elif mat_key == "metSBOTerms".casefold() or mat_key == "rxnSBOTerms".casefold():
            annotations[name] = _cell_to_str_list(
                mat_struct[caseunfold[mat_key]][0, 0], None, _sbo_re, "SBO:"
            )
        elif mat_key == "metSLMID".casefold():
            annotations[name] = _cell_to_str_list(
                mat_struct[caseunfold[mat_key]][0, 0], None, _sbo_re, "SLM:"
            )
        else:
            # If it something else, which may have commas, turn it into a list
            annotations[name] = [
                [y.strip() for y in each_cell.split(", ")] if each_cell else None
                for each_cell in _cell_to_str_list(
                    mat_struct[caseunfold[mat_key]][0, 0]
                )
            ]
    for i, obj in enumerate(target_list):
        obj.annotation = {
            prov: annotations[prov][i] for prov in providers if annotations[prov][i]
        }

    # TODO - When cobrapy.notes are revised not be a dictionary (possibly when
    #  annotations are fully SBML compliant, revise this function.


def mat_parse_notes(
    target_list: List[Object],
    mat_struct: np.ndarray,
    d_replace: str = DICT_REACTION_NOTES,
) -> None:
    """Process mat structure notes in place.

    Will process mat structured notes and add them to a list of new entities
    (metabolites, reactions, genes) in a format based on identifiers.org.

    Parameters
    ----------
    target_list: list[cobra.Object]
        A list of cobra objects, including metabolites, reactions or genes. The
        notes will be added to these lists.
    mat_struct: np.ndarray
        A darray that includes the data imported from matlab file.
    d_replace: str
        A string that points to the dictionary of converstions between MATLAB and
        notes. Default DICT_REACTION_NOTES (for reactions).
    """
    struct_names = [x.casefold() for x in mat_struct.dtype.names]
    matlab_field_dict = {
        x.casefold(): DICT_REPLACE[d_replace][x] for x in DICT_REPLACE[d_replace].keys()
    }
    caseunfold = {x.casefold(): x for x in mat_struct.dtype.names}
    annotation_matlab = list(set(struct_names).intersection(matlab_field_dict.keys()))
    note_providers = [matlab_field_dict[x] for x in annotation_matlab]
    notes = dict.fromkeys(note_providers, None)
    for name, mat_key in zip(note_providers, annotation_matlab):
        if mat_key == "rxnReferences".casefold():
            # This only removes PMID: style references. Sometimes there are other
            # things like PMC or OMIM, but those are placed as string in notes.
            _notes = _cell_to_str_list(mat_struct[caseunfold[mat_key]][0, 0])
            notes[name] = [
                _pubmed_re.sub("", x).strip()
                if x and len(_pubmed_re.sub("", x).strip())
                else None
                for x in _notes
            ]
        elif mat_key == "rxnConfidenceScores".casefold():
            notes[name] = [
                str(confidence) if confidence is not None else ""
                for confidence in _cell_to_float_list(
                    mat_struct[caseunfold[mat_key]][0, 0]
                )
            ]
        else:
            # If it something else, which may have commas, turn it into a list
            notes[name] = _cell_to_str_list(mat_struct[caseunfold[mat_key]][0, 0])
        notes[name] = [
            _punctuation_re.sub("", _double_punctuation_re.sub("", x)) if x else None
            for x in notes[name]
        ]
    for i, obj in enumerate(target_list):
        obj.notes = {prov: notes[prov][i] for prov in note_providers if notes[prov][i]}


def annotations_to_mat(
    mat_dict: OrderedDict, annotation_list: List[Dict], d_replace: str = DICT_MET_REV
) -> None:
    """Process mat structure annotations in place.

    Will process mat structured annotations and add them to a list of new entities
    (metabolites, reactions, genes) in a format based on identifiers.org.

    Parameters
    ----------
    mat_dict: OrderedDict
        An ordered dictionary having model attributes as keys and their
        respective values represented as arrays, as the values. Annotations will
        be inserted into this OrderdDict.
    annotation_list: list[Dict]
        A list of cobra annotations, in the form of a dictionary.
    d_replace: str
        A string that points to the dictionary of converstions between MATLAB and
        providers. Default DICT_MET_REV (for metabolite).
    """
    providers_used = set()
    for i in range(len(annotation_list)):
        providers_used.update(annotation_list[i].keys())
    providers_used = providers_used.intersection(DICT_REPLACE[d_replace].keys())
    providers_used = list(providers_used)
    annotation_matlab = {prov: DICT_REPLACE[d_replace][prov] for prov in providers_used}
    empty_lists = [[""] * len(annotation_list) for _ in annotation_matlab]
    annotation_cells_to_be = dict(zip(annotation_matlab.values(), empty_lists))
    for i in range(len(annotation_list)):
        for provider_key, obj_annotation in annotation_list[i].items():
            if isinstance(obj_annotation, str):
                obj_annotation = [obj_annotation]
            if provider_key == "pubmed":
                obj_annotation = ", ".join(
                    [
                        "PMID:" + annot if "PMID:" not in annot else annot
                        for annot in obj_annotation
                    ]
                )
            elif provider_key == "CHEBI":
                obj_annotation = ", ".join(
                    [
                        "CHEBI:" + annot if "CHBEI:" not in annot else annot
                        for annot in obj_annotation
                    ]
                )
            elif provider_key == "ec-code":
                obj_annotation = " or ".join(obj_annotation)
            else:
                obj_annotation = ", ".join(obj_annotation)
            if provider_key not in providers_used:
                continue
            annotation_cells_to_be[annotation_matlab[provider_key]][i] = obj_annotation
    for annotation_key, item_list in annotation_cells_to_be.items():
        mat_dict[annotation_key] = _cell(item_list)

    # TODO - When cobrapy.notes are revised not be a dictionary (possibly when
    #  annotations are fully SBML compliant, revise this function.


def notes_to_mat(
    mat_dict: OrderedDict, note_list: List[Dict], d_replace: str = DICT_MET_REV
) -> None:
    """Process mat structure notes in place.

    Will process mat structured annotations and add them to a list of new entities
    (metabolites, reactions, genes) in a format based on identifiers.org.

    Parameters
    ----------
    mat_dict: OrderedDict
        An ordered dictionary having model attributes as keys and their
        respective values represented as arrays, as the values. Annotations will
        be inserted into this OrderdDict.
    note_list: list[Dict]
        A list of cobra annotations, in the form of a dictionary.
    d_replace: str
        A string that points to the dictionary of converstions between MATLAB and
        providers. Default DICT_MET_REV (for metabolite).
    """
    providers_used = set()
    for i in range(len(note_list)):
        providers_used.update(note_list[i].keys())
    providers_used = providers_used.intersection(DICT_REPLACE[d_replace].keys())
    providers_used = list(providers_used)
    annotation_matlab = {prov: DICT_REPLACE[d_replace][prov] for prov in providers_used}
    empty_lists = [[""] * len(note_list) for _ in annotation_matlab]
    annotation_cells_to_be = dict(zip(annotation_matlab.values(), empty_lists))
    for i in range(len(note_list)):
        for provider_key, object_note in note_list[i].items():
            if provider_key not in providers_used:
                continue
            if provider_key == CONFIDENCE_STR:
                object_note = float(object_note)
            if not len(annotation_cells_to_be[annotation_matlab[provider_key]][i]):
                annotation_cells_to_be[annotation_matlab[provider_key]][i] = object_note
            else:
                # References that aren't MIRIAM compliant will go to rxnNotes
                annotation_cells_to_be[annotation_matlab[provider_key]][i] = (
                    annotation_cells_to_be[annotation_matlab[provider_key]][i]
                    + f"; {object_note}"
                )
    for annotation_key, item_list in annotation_cells_to_be.items():
        mat_dict[annotation_key] = _cell(item_list)


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
    mat["mets"] = _cell(mets.list_attr("id"))
    model_has_compartment_names = False
    if list(model.compartments.keys()) != list(model.compartments.values()):
        model_has_compartment_names = True
    if {_get_id_compartment(met_id) for met_id in mets.list_attr("id")} == {
        None
    } or model_has_compartment_names:
        comps = list(model.compartments.keys())
        mat["comps"] = _cell(comps)
        mat["compNames"] = _cell([model.compartments[comp] for comp in comps])
        mat["metComps"] = [comps.index(x) + 1 for x in mets.list_attr("compartment")]
    mat["metNames"] = _cell(mets.list_attr("name"))
    mat["metFormulas"] = _cell([str(m.formula) if m.formula else "" for m in mets])
    try:
        mat["metCharges"] = np.array(mets.list_attr("charge")) * 1.0
    except (TypeError, AttributeError):
        # can't have any None entries for charge, or this will fail
        # TODO: use custom cobra exception to handle exception
        pass
    annotations_to_mat(mat, mets.list_attr("annotation"), DICT_MET_REV)
    # TODO - When cobrapy.notes are revised not be a dictionary (possibly when
    #  annotations are fully SBML compliant, revise this function.
    notes_to_mat(mat, mets.list_attr("notes"), DICT_MET_NOTES_REV)
    mat["genes"] = _cell(model.genes.list_attr("id"))
    gene_names = model.genes.list_attr("name")
    if not all(_name == "" for _name in gene_names):
        mat["geneNames"] = _cell(gene_names)
    annotations_to_mat(mat, model.genes.list_attr("annotation"), DICT_GENE_REV)
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
    annotations_to_mat(mat, rxns.list_attr("annotation"), DICT_REACTION_REV)
    notes_to_mat(mat, rxns.list_attr("notes"), DICT_REACTION_NOTES_REV)
    rxn_subsystems = _cell(rxns.list_attr("subsystem"))
    group_ids = model.groups.list_attr("id")
    group_names = model.groups.list_attr("name")
    group_ids = [
        _id if not _name else _name for _id, _name in zip(group_ids, group_names)
    ]
    group_members_list = model.groups.list_attr("members")
    if set(rxn_subsystems) == {""} and set(group_ids) != {""}:
        subsystems = [[] for _ in rxns]
        for group_id, group_members in zip(group_ids, group_members_list):
            group = model.groups.get_by_id(group_id)
            if group.kind == "partonomy":
                for member in group_members:
                    if isinstance(member, Reaction):
                        rxn_ind = model.reactions.index(member)
                        # noinspection PyTypeChecker
                        subsystems[rxn_ind].append(group_id)
        mat["subSystems"] = _cell(
            [", ".join(subsystem_list) for subsystem_list in subsystems]
        )
    else:
        mat["subSystems"] = _cell(rxns.list_attr("subsystem"))
    stoich_mat = create_stoichiometric_matrix(model, array_type="dok")
    mat["S"] = stoich_mat if stoich_mat is not None else [[]]
    # multiply by 1 to convert to float, working around scipy bug
    # https://github.com/scipy/scipy/issues/4537
    mat["lb"] = np.array(rxns.list_attr("lower_bound")) * 1.0
    mat["ub"] = np.array(rxns.list_attr("upper_bound")) * 1.0
    mat["b"] = np.array(mets.list_attr("_bound")) * 1.0
    mat["c"] = np.array(rxns.list_attr("objective_coefficient")) * 1.0
    mat["rev"] = np.array(rxns.list_attr("reversibility")) * 1
    if model.name:
        mat["modelName"] = str(model.name)
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

    if m.dtype.names is None or {"rxns", "mets", "S", "lb", "ub"} > set(m.dtype.names):
        raise ValueError("Invalid MATLAB struct.")

    old_cobratoolbox_fields = [
        "confidenceScores",
        "metCharge",
        "ecNumbers",
        "KEGGID",
        "metSmile",
        "metHMDB",
    ]
    new_cobratoolbox_fields = [
        "rxnConfidenceScores",
        "metCharges",
        "rxnECNumbers",
        "metKEGGID",
        "metSmiles",
        "metHMDBID",
    ]
    for old_field, new_field in zip(old_cobratoolbox_fields, new_cobratoolbox_fields):
        if old_field in m.dtype.names and new_field not in m.dtype.names:
            logger.warning(
                f"This model seems to have {old_field} instead of {new_field} field. "
                f"Will use {old_field} for what {new_field} represents."
            )
            new_names = list(m.dtype.names)
            new_names[new_names.index(old_field)] = new_field
            m.dtype.names = new_names

    model = Model()
    if model_id is not None:
        model.id = model_id
    elif "description" in m.dtype.names:
        description = m["description"][0, 0][0]
        if not isinstance(description, str) and len(description) > 1:
            model.id = description[0]
            logger.warning("Several IDs detected, only using the first.")
        else:
            model.id = description
    else:
        model.id = "imported_model"
    if "modelName" in m.dtype.names and np.size(m["modelName"][0, 0]):
        model.name = m["modelName"][0, 0][0]

    met_ids = _cell_to_str_list(m["mets"][0, 0])
    if {"metComps", "comps", "compNames"}.issubset(m.dtype.names):
        met_comp_index = [x[0] - 1 for x in m["metComps"][0][0]]
        comps = _cell_to_str_list(m["comps"][0, 0])
        comp_names = _cell_to_str_list(m["compNames"][0][0])
        met_comps = [comps[i] for i in met_comp_index]
        met_comp_names = [comp_names[i] for i in met_comp_index]
    else:
        logger.warning(
            f"No defined compartments in model {model.id}. "
            f"Compartments will be deduced heuristically "
            f"using regular expressions."
        )
        met_comps = [_get_id_compartment(x) for x in met_ids]
        met_comp_names = met_comps
        if None in met_comps or "" in met_comps:
            raise ValueError("Some compartments were empty. Check the model!")
        logger.warning(
            f"Using regular expression found the following compartments:"
            f"{', '.join(sorted(set(met_comps)))}"
        )
    if None in met_comps or "" in met_comps:
        raise ValueError("Some compartments were empty. Check the model!")
    model.compartments = dict(zip(met_comps, met_comp_names))
    met_names, met_formulas, met_charges = None, None, None
    try:
        met_names = _cell_to_str_list(m["metNames"][0, 0], "")
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
    new_metabolites = []
    for i in range(len(met_ids)):
        new_metabolite = Metabolite(met_ids[i], compartment=met_comps[i])
        if met_names:
            new_metabolite.name = met_names[i]
        if met_charges:
            new_metabolite.charge = met_charges[i]
        if met_formulas:
            new_metabolite.formula = met_formulas[i]
        new_metabolites.append(new_metabolite)
    mat_parse_annotations(new_metabolites, m, d_replace=DICT_MET)
    mat_parse_notes(new_metabolites, m, d_replace=DICT_MET_NOTES)
    model.add_metabolites(new_metabolites)

    if "genes" in m.dtype.names:
        gene_names = None
        gene_ids = _cell_to_str_list(m["genes"][0, 0])
        try:
            gene_names = _cell_to_str_list(m["geneNames"][0, 0])
        except (IndexError, ValueError):
            # TODO: use custom cobra exception to handle exception
            pass
        new_genes = [
            Gene(gene_ids[i], name=gene_names[i]) if gene_names else Gene(gene_ids[i])
            for i in range(len(gene_ids))
        ]
        mat_parse_annotations(new_genes, m, d_replace=DICT_GENE)
        for current_gene in new_genes:
            current_gene._model = model
        model.genes += new_genes

    new_reactions = []
    rxn_ids = _cell_to_str_list(m["rxns"][0, 0])
    rxn_lbs = _cell_to_float_list(m["lb"][0, 0], empty_value=None, inf_value=inf)
    rxn_ubs = _cell_to_float_list(m["ub"][0, 0], empty_value=None, inf_value=inf)
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
        # RECON3.0 mat has an array within an array for subsystems.
        # If we find a model that has multiple subsytems per reaction, this should be
        # modified
        if np.sctype2char(m["subSystems"][0, 0][0][0]) == "O" and isinstance(
            m["subSystems"][0, 0][0][0][0], np.ndarray
        ):
            rxn_subsystems = [
                each_cell[0][0][0][0] if each_cell else ""
                for each_cell in m["subSystems"][0, 0]
            ]
        # Other matlab files seem normal.
        else:
            rxn_subsystems = _cell_to_str_list(m["subSystems"][0, 0], "")
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
    mat_parse_annotations(new_reactions, m, d_replace=DICT_REACTION)
    # TODO - When cobrapy.notes are revised not be a dictionary (possibly when
    #  annotations are fully SBML compliant, revise this function.
    mat_parse_notes(new_reactions, m, d_replace=DICT_REACTION_NOTES)

    csc = scipy_sparse.csc_matrix(m["S"][0, 0])
    for i in range(csc.shape[1]):
        stoic_dict = {
            model.metabolites[j]: csc[j, i] for j in csc.getcol(i).nonzero()[0]
        }
        new_reactions[i].add_metabolites(stoic_dict)

    model.add_reactions(new_reactions)

    # Make subsystems into groups
    if rxn_subsystems:
        rxn_group_names = set(rxn_subsystems).difference({None})
        new_groups = []
        for g_name in sorted(rxn_group_names):
            group_members = model.reactions.query(
                lambda x: x.subsystem == g_name  # noqa: B023
            )
            new_group = Group(
                id=g_name, name=g_name, members=group_members, kind="partonomy"
            )
            new_group.annotation["sbo"] = "SBO:0000633"
            new_groups.append(new_group)
        model.add_groups(new_groups)

    if "c" in m.dtype.names:
        c_vec = _cell_to_float_list(m["c"][0, 0])
        coefficients = dict(zip(new_reactions, c_vec))
        set_objective(model, coefficients)
    else:
        logger.warning("Objective vector `c` not found.")

    if "osenseStr" in m.dtype.names:
        if isinstance(m["osenseStr"][0, 0][0], np.str_):
            model.objective_direction = str(m["osenseStr"][0, 0][0])
        elif isinstance(m["osenseStr"][0, 0][0], np.ndarray):
            model.objective_direction = str(m["osenseStr"][0, 0][0][0])
    elif "osense" in m.dtype.names:
        osense = float(m["osense"][0, 0][0][0])
        objective_direction_str = "max"
        if osense == 1:
            objective_direction_str = "min"
        model.objective_direction = objective_direction_str
    return model
