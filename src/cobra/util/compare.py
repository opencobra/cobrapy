"""Comparing models, reactions, metabolites, genes and groups."""

from typing import TYPE_CHECKING, Dict, Optional, Tuple, TypeVar


if TYPE_CHECKING:
    from cobra import Group, Model, Object, Reaction

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
    same = set(o for o in shared_keys if d1[o] == d2[o])
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

    Parameters
    ----------
    model1: cobra.Model
        Model to compare.
    model2: cobra.Model
        Other Model to compare.
    ignore_notes: bool, optional
        Whether or not to ignore the notes field in the
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
    model_comparison["metabolites"] = dict()
    model_comparison["different_mets"] = list()
    mets_model1 = set(model1.metabolites.list_attr("id"))
    mets_model2 = set(model2.metabolites.list_attr("id"))
    if mets_model1 != mets_model2:
        if mets_model1 - mets_model2:
            model_comparison["metabolites"]["added"] = mets_model1 - mets_model2
        if mets_model2 - mets_model1:
            model_comparison["metabolites"]["removed"] = mets_model2 - mets_model1
    for _id in list(mets_model1.intersection(mets_model2)):
        met1 = model1.metabolites.get_by_id(_id)
        met2 = model2.metabolites.get_by_id(_id)
        _eq, _comparison = compare_state(met1, met2, ignore_keys=ignore_keys)
        if not _eq:
            _is_equivalent = False
            model_comparison["different_mets"].append(_id)
        model_comparison["metabolites"][_id] = _comparison

    model_comparison["reactions"] = dict()
    model_comparison["different_rxns"] = list()
    rxns_model1 = set(model1.reactions.list_attr("id"))
    rxns_model2 = set(model2.reactions.list_attr("id"))
    if rxns_model1 - rxns_model2:
        model_comparison["reactions"]["added"] = rxns_model1 - rxns_model2
    if rxns_model2 - rxns_model1:
        model_comparison["reactions"]["removed"] = rxns_model2 - rxns_model1
    for _id in list(rxns_model1.intersection(rxns_model2)):
        rxn1 = model1.reactions.get_by_id(_id)
        rxn2 = model2.reactions.get_by_id(_id)
        _eq, _comparison = compare_reaction_state(rxn1, rxn2, ignore_keys=ignore_keys)
        if not _eq:
            _is_equivalent = False
            model_comparison["different_rxns"].append(_id)
        model_comparison["reactions"][_id] = _comparison

    model_comparison["genes"] = dict()
    model_comparison["different_genes"] = list()
    genes_model1 = set(model1.genes.list_attr("id"))
    genes_model2 = set(model2.genes.list_attr("id"))
    if genes_model1 - genes_model2:
        model_comparison["genes"]["added"] = genes_model1 - genes_model2
    if genes_model2 - genes_model1:
        model_comparison["genes"]["removed"] = genes_model2 - genes_model1
    for _id in list(genes_model1.intersection(genes_model2)):
        gene1 = model1.genes.get_by_id(_id)
        gene2 = model2.genes.get_by_id(_id)
        _eq, _comparison = compare_state(gene1, gene2, ignore_keys=ignore_keys)
        if not _eq:
            _is_equivalent = False
            model_comparison["different_genes"].append(_id)
        model_comparison["genes"][_id] = _comparison

    model_comparison["groups"] = dict()
    model_comparison["different_groups"] = list()
    groups_model1 = set(model1.groups.list_attr("id"))
    groups_model2 = set(model2.groups.list_attr("id"))
    if groups_model1 - groups_model2:
        model_comparison["groups"]["added"] = groups_model1 - groups_model2
    if groups_model2 - groups_model1:
        model_comparison["groups"]["removed"] = groups_model2 - groups_model1
    for _id in list(groups_model1.intersection(groups_model2)):
        group1 = model1.groups.get_by_id(_id)
        group2 = model2.groups.get_by_id(_id)
        _eq, _comparison = compare_state(group1, group2, ignore_keys=ignore_keys)
        if not _eq:
            _is_equivalent = False
            model_comparison["different_groups"].append(_id)
        model_comparison["groups"][_id] = _comparison

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
