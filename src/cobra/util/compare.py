"""Comparing models, reactions, metabolites, genes and groups."""

from typing import (
    TYPE_CHECKING,
    Callable,
    Dict,
    List,
    Optional,
    Set,
    Tuple,
    TypeVar,
    Union,
)


if TYPE_CHECKING:
    from cobra.core import DictList, Group, Model, Object, Reaction

    TObject = TypeVar("TObject", bound=Object)


def dict_compare(
    d1: Dict, d2: Dict, _dont_compare: Optional[set] = None
) -> Dict[str, Union[Set, Dict]]:
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

    Returns
    -------
    A dicitonary comprised of
        "added" - set of attributes present in the first, but not the second
        "removed" - set of attributes present in the second, but not the first
        "same" - dict of keys prsent in both with the same values
        "modified" - dict of keys present in both with different values. Each key will
                    have a tuple of values in both dictionaries given.
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
    """Cmpare two cobra Objects (including subclasses).

    Useful for Metaboite, and Gene Comparison.
    Not useful for comparing GPRs(). Use the equality in GPRs() directly.
    For Reaction and Group, use the specific functions which do some pre-processing.

    Parameters
    ----------
    obj1: Object, Metabolite, Gene
    obj2: Object, Metabolite, Gene
    ignore_keys: Set, optional
        A set of keys to ignore. Defuault None (empty set - all keys will be compared).

    Returns
    -------
    Tuple - bool, Dict
        A tuple of a boolean (are the two objects equivalent or different) and
        a dictionary specifying how they differed.
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
    """Compare two cobra Reactions.

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
        A tuple of a boolean (are the two objects equivalent or different) and a
        dictionary specifying how they differed.
    """
    if ignore_keys is None:
        ignore_keys = set()
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
    """Compare two cobra Groups.

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
        A tuple of a boolean (are the two objects equivalent or different) and a
        dictionary specifying how they differed.
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
    dictlist1: "DictList",
    dictlist2: "DictList",
    ignore_keys: Optional[set] = None,
    comparison_function: Optional[Callable] = None,
) -> Tuple[bool, Dict, List]:
    """Compare dictlist of objects. Useful when comparing models.

    Will check whether there are objects in dictlist1 that aren't present in dictlist2,
    and vice versa.
    Objects that appear in both dictlists will be compared using the comparison
    function, which allows different functions for reaction or group. The default is
    None, which will result in simple compare_state() used as the comparison function.

    Parameters
    ----------
    ignore_keys: set, optional
        What keys should be ignored. Defualt None.
    dictlist1: cobra.DictList
        The first dictlist to compare
    dictlist2: cobra.DictList
    comparison_function: Callable
        A function to use for comparing the objects in the dictlists.

    Returns
    -------
    Tuple: bool, Dict, List
        A tuple of a boolean (are the two dictlists equivalent or different) and a
        dictionary specifying how they differed. It also returns a list of ids that
        were present in both dictlists, but different.

    See Also
    --------
    compare_state()
    compare_reaction_state()
    compare_group_state()
    """
    if comparison_function is None:
        comparison_function = compare_state
    _is_equivalent = True
    comparison = {}
    diff_objs = []
    ids_model1 = set(dictlist1.list_attr("id"))
    ids_model2 = set(dictlist2.list_attr("id"))
    if ids_model1 - ids_model2:
        comparison["added"] = ids_model1 - ids_model2
        _is_equivalent = False
    if ids_model2 - ids_model1:
        comparison["removed"] = ids_model2 - ids_model1
        _is_equivalent = False
    for _id in dictlist1.intersection(dictlist2):
        obj1 = dictlist1.get_by_id(_id)
        obj2 = dictlist2.get_by_id(_id)
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
) -> Tuple[bool, Dict]:
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
    ignore_keys: set, optional
        What keys should be ignored. Defualt None.

    Returns
    -------
    Tuple - bool, Dict
        A tuple of a boolean (are the two models equivalent or different)
        and a dictionary specifying how they differed. The dictionary contains
        different_x as a list and x for metabolites, reactions, genes, groups.
        The differenet_x specifies which comparisons were not equivalent, while the
        x contains the full dictionary of comparing each element (each group,
        metabolite, reaction, gene).

    See Also
    --------
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
        model1.metabolites, model2.metabolites, ignore_keys
    )
    model_comparison["metabolites"] = comp_result
    model_comparison["different_mets"] = different_ids
    _is_equivalent &= _eq

    _eq, comp_result, different_ids = compare_dictlists(
        model1.reactions,
        model2.reactions,
        ignore_keys,
        comparison_function=compare_reaction_state,
    )
    model_comparison["reactions"] = comp_result
    model_comparison["different_reactions"] = different_ids
    _is_equivalent &= _eq

    _eq, comp_result, different_ids = compare_dictlists(
        model1.genes, model2.genes, ignore_keys
    )
    model_comparison["genes"] = comp_result
    model_comparison["different_genes"] = different_ids
    _is_equivalent &= _eq

    _eq, comp_result, different_ids = compare_dictlists(
        model1.groups,
        model2.groups,
        ignore_keys,
        comparison_function=compare_group_state,
    )
    model_comparison["genes"] = comp_result
    model_comparison["different_genes"] = different_ids
    _is_equivalent &= _eq

    return _is_equivalent, model_comparison


def fix_for_reaction_notes_changes(diff_dict: Dict, diff_list: List) -> None:
    r"""Fix differences caused in reaction Notes when reading and writing models.

    If you wish to compare reaction Notes, there may be some changes that are caused
    because the reading and writing functions (in Matlab) do not fully preserve all
    the info in the Notes field.

    If you're using this on a model comparison, it should be the two fields
    'reactions' and 'different_reactions'.

    Matlab -    reaction Notes may include References as well. These are read as Notes,
                outputted to the rxnReferences field if they match the pubmed format
                r"PMID: ?\d+"

    Parameters
    ----------
    diff_dict: Dict
        Dictionary of differences.
    diff_list: List
        List of reactions that were present in both dictlists/models but had were
        different in some fields.

    Examples
    --------
    >>>>
    >>>> comparison = compare_model_state(model1, model2)

    """
    for key in diff_list:
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
                diff_list.remove(key)
                diff_dict[key]["modified"].__delitem__("notes")
