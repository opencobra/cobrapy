"""Comparing models, reactions, metabolites, genes and groups."""

from typing import TYPE_CHECKING, Dict, Optional, Tuple, TypeVar


if TYPE_CHECKING:
    from cobra.core import Group, Model, Object, Reaction

    TObject = TypeVar("TObject", bound=Object)


def dict_compare(d1: Dict, d2: Dict, _dont_compare: Optional[set] = None):
    """Compare two dictionaries.

    This function will identify overlapping keys, added, removed keys between
    dictonaries. If there are identical keys which will not have the same value, they
    will be noted as 'modified'.

    Parameters
    ----------
    d1: dict
        Dictionary to compare.
    d2: dict
        Dictionary to compare.
    _dont_compare: set
        Keys that should not be compared. Optional. Default None (compare all keys).
    """
    if _dont_compare is None:
        _dont_compare = set()
    d1_keys = set(d1.keys()).difference(_dont_compare)
    d2_keys = set(d2.keys()).difference(_dont_compare)
    shared_keys = d1_keys.intersection(d2_keys)
    added = d1_keys - d2_keys
    removed = d2_keys - d1_keys
    modified = {o: (d1[o], d2[o]) for o in shared_keys if d1[o] != d2[o]}
    same = {o for o in shared_keys if d1[o] == d2[o]}
    return {"added": added, "removed": removed, "modified": modified, "same": same}


def compare_state(
    obj1: "TObject", obj2: "TObject", ignore_keys: Optional[set] = None
) -> Tuple[bool, Dict]:
    """Will compare two cobra Objects (and what is derived from them).

    Not useful for comparing GPRs(). Use the equality in GPRs() directly.
    For Reaction and Group, use the specific functions which do some processing.

    Parameters
    ----------
    obj1: Object, Metabolite, Gene
    obj2: Object, Metabolite, Gene
    ignore_keys: Set, optional
        A set of keys to ignore. Defuault None (empty set - all keys will be compared).

    Returns
    -------
    Tuple - bool, Dict
        A tuple of a boolean (are the two objects different or not) and a dictionary
        specifying how they differed.
    """
    _is_equivalent = True
    if ignore_keys is None:
        ignore_keys = set()
    _is_equivalent = True
    state1 = obj1.__getstate__()
    state2 = obj2.__getstate__()
    _comparison = dict_compare(state1, state2, ignore_keys)
    if _comparison["added"] or _comparison["removed"] or _comparison["modified"]:
        _is_equivalent = False
    return _is_equivalent, _comparison


def compare_reaction_state(
    rxn1: "Reaction", rxn2: "Reaction", ignore_keys: Optional[set] = None
) -> Tuple[bool, Dict]:
    """Will compare two cobra Reactions.

    In order to avoid recursion and disagreement on memory address
    genes are transformed to gene.ids
    metabolites are transformed to metabolite.ids

    Parameters
    ----------
    rxn1: Reaction
    rxn2: Reaction
    ignore_keys: Set, optional
        A set of keys to ignore. Defuault None (empty set - all keys will be compared).

    Returns
    -------
    Tuple - bool, Dict
        A tuple of a boolean (are the two objects different or not) and a dictionary
        specifying how they differed.
    """
    _is_equivalent = True
    state1 = rxn1.__getstate__()
    state1["_metabolites"] = {met.id: stoic for met, stoic in rxn1.metabolites.items()}
    state1["_genes"] = {gene.id for gene in rxn1.genes}
    state2 = rxn2.__getstate__()
    state2["_metabolites"] = {met.id: stoic for met, stoic in rxn2.metabolites.items()}
    state2["_genes"] = {gene.id for gene in rxn2.genes}
    _comparison = dict_compare(state1, state2, ignore_keys)
    if _comparison["added"] or _comparison["removed"] or _comparison["modified"]:
        _is_equivalent = False
    return _is_equivalent, _comparison


def compare_group_state(
    group1: "Group", group2: "Group", ignore_keys: Optional[set] = None
) -> Tuple[bool, Dict]:
    """Will compare two cobra Groups.

    Members are transformed to a list of reaction ids in order to avoid differences in
    memory address leading to false positives.

    Parameters
    ----------
    group1: Group
    group2: Group
    ignore_keys: Set, optional
        A set of keys to ignore. Defuault None (empty set - all keys will be compared).

    Returns
    -------
    Tuple - bool, Dict
        A tuple of a boolean (are the two objects different or not) and a dictionary
        specifying how they differed.
    """
    _is_equivalent = True
    state1 = group1.__getstate__()
    state2 = group2.__getstate__()
    state1["_members"] = group1.members.list_attr("id")
    state2["_members"] = group2.members.list_attr("id")
    _comparison = dict_compare(state1, state2, ignore_keys)
    if _comparison["added"] or _comparison["removed"] or _comparison["modified"]:
        _is_equivalent = False
    return _is_equivalent, _comparison


def compare_dictlists(
    ignore_keys, dictlist_model1, dict_list_model2, comparison_function=None
):
    if comparison_function is None:
        comparison_function = compare_state
    _is_equivalent = True
    comparison = {}
    diff_objs = []
    ids_model1 = set(dictlist_model1.list_attr("id"))
    ids_model2 = set(dict_list_model2.list_attr("id"))
    if ids_model1 - ids_model2:
        comparison["added"] = ids_model1 - ids_model2
        _is_equivalent = False
    if ids_model2 - ids_model1:
        comparison["removed"] = ids_model2 - ids_model1
        _is_equivalent = False
    for _id in dictlist_model1.intersection(dict_list_model2):
        obj1 = dictlist_model1.get_by_id(_id)
        obj2 = dict_list_model2.get_by_id(_id)
        _eq, _comparison = comparison_function(obj1, obj2, ignore_keys=ignore_keys)
        if not _eq:
            _is_equivalent = False
            diff_objs.append(_id)
        comparison[_id] = _comparison
    return _is_equivalent, comparison, diff_objs


def compare_model_state(
    model1: "Model",
    model2: "Model",
    ignore_notes: bool = True,
    ignore_keys: Optional[set] = None,
):
    """Recursively compare model states.

    Will compare the model and then compare metabolites, reactions, genes, groups in
    the model. Models will be considered different if any of the objects within the
    cobra model are different.
    Will ignore the notes field. If you want to compare the notes field, you may need
    to run _fix_xml_annotation_to_identifiers and/or fix_for_notes_changes.

    Parameters
    ----------
    model1: cobra.Model
        Model to compare.
    model2: cobra.Model
        Other Model to compare.
    ignore_notes: bool, optional
        Whether or not to ignore the notes field in the model. Default True.
    ignore_keys

    Returns
    -------
    Tuple - bool, Dict
        A tuple of a boolean (are the two models different or not)
        and a dictionary specifying how they differed. The dictionary contains
        different_x as a list and x for metabolites, reactions, genes, groups.
        The differenet_x specifies which comparisons were not equivalent, while the
        x contains the full dictionary of comparing each element (each group,
        metabolite, reaction, gene).

    See Also
    --------
    _fix_xml_annotation_to_identifiers()
    fix_for_notes_changes()
    """
    _is_equivalent = True
    if ignore_keys is None:
        ignore_keys = set()
    if ignore_notes:
        ignore_keys = ignore_keys.union({"notes"})
    do_not_compare_models = {
        "metabolites",
        "reactions",
        "genes",
        "notes",
        "annotation",
        "_annotation",
        "groups",
        "_sbml",  # don't care about SBML properties of the file, just how it is read
        "_id",  # Will often be different based on different files
        "_solver",  # Will be different memory locations
    }
    _eq, model_comparison = compare_state(model1, model2, do_not_compare_models)
    _is_equivalent = _eq
    _eq, comp_result, different_ids = compare_dictlists(
        ignore_keys, model1.metabolites, model2.metabolites
    )
    model_comparison["metabolites"] = comp_result
    model_comparison["different_mets"] = different_ids
    _is_equivalent &= _eq

    _eq, comp_result, different_ids = compare_dictlists(
        ignore_keys,
        model1.reactions,
        model2.reactions,
        comparison_function=compare_reaction_state,
    )
    model_comparison["reactions"] = comp_result
    model_comparison["different_reactions"] = different_ids
    _is_equivalent &= _eq

    _eq, comp_result, different_ids = compare_dictlists(
        ignore_keys, model1.genes, model2.genes
    )
    model_comparison["genes"] = comp_result
    model_comparison["different_genes"] = different_ids
    _is_equivalent &= _eq

    _eq, comp_result, different_ids = compare_dictlists(
        ignore_keys,
        model1.groups,
        model2.groups,
        comparison_function=compare_group_state,
    )
    model_comparison["genes"] = comp_result
    model_comparison["different_genes"] = different_ids
    _is_equivalent &= _eq

    return _is_equivalent, model_comparison


def _fix_xml_annotation_to_identifiers(model: "Model") -> None:
    """Fix XML models which have annotations that do not match identifiers.org.

    This function will fix the dict keys of annotations to match identifiers.org.
    Eventually, the XML models should be fixed and cobrapy should be strict, but this is
    part of SBML rewriting of annotations
    see: https://github.com/opencobra/cobrapy/issues/684

    Useful for comapring matlab models with XML models, otherwise the difference in
    annotation behavoir confuses the funciton.

    Parameters
    ----------
    model: Model
        A model to fix
    """
    for met in model.metabolites:
        if met.formula == "":
            met.formula = None
        if len(met.annotation):
            if "chebi" in met.annotation.keys():
                met.annotation["CHEBI"] = met.annotation.pop("chebi")
            if "sbo" in met.annotation.keys():
                met.annotation["SBO"] = met.annotation.pop("sbo")
            for annot, val in met.annotation.items():
                if isinstance(val, str):
                    met.annotation[annot] = [val]
    for rxn in model.reactions:
        rxn.name = rxn.name.strip()
        if "sbo" in rxn.annotation.keys():
            rxn.annotation["SBO"] = rxn.annotation.pop("sbo")
        if len(rxn.annotation):
            for annot, val in rxn.annotation.items():
                if isinstance(val, str):
                    rxn.annotation[annot] = [val]
    for gene in model.genes:
        if len(gene.annotation):
            if "ncbigi" in gene.annotation.keys():
                gene.annotation["ncbiprotein"] = gene.annotation.pop("ncbigi")
            for annot, val in gene.annotation.items():
                if isinstance(val, str):
                    gene.annotation[annot] = [val]


def fix_for_notes_changes(diff_dict, diff_set):
    for key in list(diff_set):
        if "notes" in diff_dict[key]["modified"].keys():
            note_dictionaries = diff_dict[key]["modified"]["notes"]
            note_dictionaries[0] = {
                k: v
                for k, v in note_dictionaries[0].items()
                if k != "References" and k != "NOTES"
            }
            note_dictionaries[1] = {
                k: v
                for k, v in note_dictionaries[1].items()
                if k != "References" and k != "NOTES"
            }
            if note_dictionaries[0] == note_dictionaries[1]:
                diff_set.remove(key)
                diff_dict[key]["modified"].__delitem__("notes")
