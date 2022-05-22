"""Provide functions for cobrapy objects to generic Python objects and vice-versa."""
import itertools
import re
from collections import OrderedDict, defaultdict
from typing import TYPE_CHECKING, Dict, List, Sequence, Set, Tuple, Union

import numpy as np

from ..core import Gene, Group, Metabolite, Model, Reaction
from ..core.metadata import MetaData, Notes
from ..core.metadata.helper import (
    URL_IDENTIFIERS_PATTERN,
    parse_identifiers_uri,
)
from ..io.sbml import (
    F_GENE,
    F_GENE_REV,
    F_GROUP,
    F_GROUP_REV,
    F_REACTION,
    F_REACTION_REV,
    F_REPLACE,
    F_SPECIE,
    F_SPECIE_REV,
)
from ..util.solver import set_objective


if TYPE_CHECKING:
    from cobra import Object

_REQUIRED_REACTION_ATTRIBUTES = [
    "id",
    "name",
    "lower_bound",
    "upper_bound",
    "gene_reaction_rule",
]
_ORDERED_OPTIONAL_REACTION_KEYS = [
    "objective_coefficient",
    "subsystem",
    "notes",
    "annotation",
]
_OPTIONAL_REACTION_ATTRIBUTES = {
    "objective_coefficient": 0,
    "subsystem": "",
    "notes": {},
    "annotation": {},
}

_REQUIRED_METABOLITE_ATTRIBUTES = ["id", "name", "compartment"]
_ORDERED_OPTIONAL_METABOLITE_KEYS = [
    "charge",
    "formula",
    "_bound",
    "notes",
    "annotation",
]
_OPTIONAL_METABOLITE_ATTRIBUTES = {
    "charge": None,
    "formula": None,
    "_bound": 0,
    "notes": {},
    "annotation": {},
}

_REQUIRED_GENE_ATTRIBUTES = ["id", "name"]
_ORDERED_OPTIONAL_GENE_KEYS = ["notes", "annotation"]
_OPTIONAL_GENE_ATTRIBUTES = {
    "notes": {},
    "annotation": {},
}

_REQUIRED_GROUP_ATTRIBUTES = ["id", "kind", "members"]
_ORDERED_OPTIONAL_GROUP_KEYS = ["name", "notes", "annotation"]
_OPTIONAL_GROUP_ATTRIBUTES = {
    "name": "",
    "notes": {},
    "annotation": {},
}

_ORDERED_OPTIONAL_MODEL_KEYS = ["name", "compartments", "notes", "annotation"]
_OPTIONAL_MODEL_ATTRIBUTES = {
    "name": None,
    #  "description": None, should not actually be included
    "compartments": [],
    "notes": {},
    "annotation": {},
}


def flatten(list_of_lists: Union[List, Tuple]) -> List:
    """Flatten lists of lists.

    Parameters
    ----------
    list_of_lists: List or Tuple
        List or Tuple of lists or tuples to flatten.

    Returns
    -------
    list: flattened list

    """
    return list(itertools.chain.from_iterable(list_of_lists))


def _fix_type(
    value: Union[str, np.float, np.bool, Set, Dict]
) -> Union[str, float, bool, List, OrderedDict, Dict]:
    """Convert possible types to correct Python types.

    Parameters
    ----------
    value : str, np.float, np.bool, set, dict
        The value to fix type for.

    Returns
    -------
    str, float, bool, list, dict
        The fixed type for the value.

    """
    # Because numpy floats can not be pickled to json
    if isinstance(value, str):
        return str(value)
    if isinstance(value, np.float) and (np.isnan(value) or np.isinf(value)):
        return str(value)
    if isinstance(value, np.float):
        return float(value)
    if isinstance(value, np.bool):
        return bool(value)
    if isinstance(value, set):
        return list(value)
    if isinstance(value, dict):
        return OrderedDict((key, value[key]) for key in sorted(value))
    if isinstance(value, Notes):
        return str(value)
    if isinstance(value, MetaData):
        return value.to_dict()
    # handle legacy Formula type
    if value.__class__.__name__ == "Formula":
        return str(value)
    if value is None:
        return ""
    return value


def _update_optional(
    cobra_object: "Object",
    new_dict: Dict,
    optional_attribute_dict: Dict,
    ordered_keys: Sequence,
) -> None:
    """Update `new_dict` with optional attributes from `cobra_object`.

    Parameters
    ----------
    cobra_object : cobra.Object
        The cobra Object to update optional attributes from.
    new_dict : dict
        The dictionary to update optional attributes for.
    optional_attribute_dict : dict
        The dictionary to use as default value lookup store.
    ordered_keys : list, tuple
        The keys to update values for.

    Raises
    ------
    IndexError
        If key in `ordered_keys` is not found in `optional_attribute_dict`.
    AttributeError
        If key in `ordered_keys` is not found in `cobra_object`.

    """
    for key in ordered_keys:
        default = optional_attribute_dict[key]
        value = getattr(cobra_object, key)
        if key == "notes" and (
            value is None or value.notes_xhtml is None or len(value.notes_xhtml) == 0
        ):
            continue
        if value is None or value == default:
            continue
        new_dict[key] = _fix_type(value)


def _fix_id_from_dict(
    _id_to_fix: str,
    _class_to_fix_to: str,
    f_replace: dict = F_REPLACE,  # noqa:    W0102
):
    if f_replace is None:
        f_replace = {}
    if not f_replace:
        return _id_to_fix
    if _class_to_fix_to == "Metabolite":
        return F_REPLACE[F_SPECIE](_id_to_fix)
    elif _class_to_fix_to == "Reaction":
        return F_REPLACE[F_REACTION](_id_to_fix)
    elif _class_to_fix_to == "Gene":
        return F_REPLACE[F_GENE](_id_to_fix)
    elif _class_to_fix_to == "Group":
        return F_REPLACE[F_GROUP](_id_to_fix)


def _fix_value_from_dict(_key: str, _value_to_fix: Union[List, str]):
    if _key == "annotation":
        # if annotation is in the form of list of list, modify the format
        # https://github.com/opencobra/cobrapy/issues/736
        if isinstance(_value_to_fix, list) and isinstance(_value_to_fix[0], list):
            _value_to_fix = flatten(_value_to_fix)
        anno_dict = defaultdict(list)
        if isinstance(_value_to_fix, list):
            for item in _value_to_fix:
                if re.match(URL_IDENTIFIERS_PATTERN, item):
                    provider, identifier = parse_identifiers_uri(item)
                    anno_dict[provider].append(identifier)
            _value_to_fix = anno_dict
        _value_to_fix = MetaData.from_dict(_value_to_fix)
    elif _key == "notes":
        _value_to_fix = Notes(_value_to_fix)
    elif _key == "lower_bound" or _key == "upper_bound":
        _value_to_fix = float(_value_to_fix)

    return _value_to_fix


def _get_by_id(
    _id: str, _object_to_get: str, model: Model
) -> Union[Gene, Metabolite, Group, Reaction]:
    if _object_to_get == "Reaction":
        return model.reactions.get_by_id(_id)
    elif _object_to_get == "Metabolite":
        return model.metabolites.get_by_id(_id)
    elif _object_to_get == "Group":
        return model.groups.get_by_id(_id)
    elif _object_to_get == "Gene":
        return model.genes.get_by_id(_id)


def _metabolite_to_dict(
    metabolite: Metabolite, f_replace: dict = F_REPLACE  # noqa:    W0102
) -> OrderedDict:
    """Convert a cobra Metabolite object to dictionary.

    Parameters
    ----------
    metabolite : cobra.Metabolite
        The cobra.Metabolite to convert to dictionary.
    f_replace : dict of replacement functions for id replacement
        Dictionary of replacement functions for gene, specie, and reaction.
        By default, the following id changes are performed on import:
        clip G_ from genes, clip M_ from species, clip R_ from reactions
        If no replacements should be performed, set f_replace={} or None

    Returns
    -------
    dict
        The converted dictionary object.

    See Also
    --------
    _metabolite_from_dict : Convert a dictionary to cobra Metabolite object.

    """
    if f_replace is None:
        f_replace = {}

    new_metabolite = OrderedDict()
    for key in _REQUIRED_METABOLITE_ATTRIBUTES:
        if key == "id" and f_replace and F_SPECIE_REV in f_replace:
            new_metabolite[key] = _fix_type(f_replace[F_SPECIE_REV](metabolite.id))
        else:
            new_metabolite[key] = _fix_type(getattr(metabolite, key))
    _update_optional(
        metabolite,
        new_metabolite,
        _OPTIONAL_METABOLITE_ATTRIBUTES,
        _ORDERED_OPTIONAL_METABOLITE_KEYS,
    )
    return new_metabolite


def _metabolite_from_dict(
    metabolite: Dict, f_replace: dict = F_REPLACE # noqa:    W0102
) -> Metabolite:
    """Convert a dictionary to cobra Metabolite object.

    Parameters
    ----------
    metabolite : dict
        The dictionary to convert to cobra.Metabolite .

    Returns
    -------
    cobra.Metabolite

    See Also
    --------
    _metabolite_to_dict : Convert a cobra Metabolite object to dictionary.

    """
    if f_replace is None:
        f_replace = {}

    new_metabolite = Metabolite(
        _fix_id_from_dict(metabolite["id"], "Metabolite", f_replace)
    )
    [
        setattr(new_metabolite, k, _fix_value_from_dict(k, v))
        for k, v in metabolite.items()
        if k != "id"
    ]
    return new_metabolite


def _gene_to_dict(gene: Gene) -> OrderedDict:
    """Convert a cobra Gene object to dictionary.

    Parameters
    ----------
    gene : cobra.Gene
        The cobra.Gene to convert to dictionary.

    Returns
    -------
    dict
        The converted dictionary object.

    See Also
    --------
    _gene_from_dict : Convert a dictionary to cobra Gene object.

    """
    new_gene = OrderedDict()
    for key in _REQUIRED_GENE_ATTRIBUTES:
        if key == "id":
            new_gene[key] = _fix_type(F_REPLACE[F_GENE_REV](gene.id))
        else:
            new_gene[key] = _fix_type(getattr(gene, key))
    _update_optional(
        gene, new_gene, _OPTIONAL_GENE_ATTRIBUTES, _ORDERED_OPTIONAL_GENE_KEYS
    )
    return new_gene


def gene_from_dict(gene: Dict, f_replace: dict = F_REPLACE) -> Gene:  # noqa:    W0102
    """Convert a dictionary to cobra Gene object.

    Parameters
    ----------
    gene : dict
        The dictionary to convert to cobra.Gene .

    Returns
    -------
    cobra.Gene
        The converted cobra.Gene object.

    See Also
    --------
    _gene_to_dict : Convert a cobra Gene object to a dictionary.

    """
    if f_replace is None:
        f_replace = {}

    new_gene = Gene(_fix_id_from_dict(gene["id"], "Gene", f_replace))
    [
        setattr(new_gene, k, _fix_value_from_dict(k, v))
        for k, v in gene.items()
        if k != "id"
    ]
    return new_gene


def _reaction_to_dict(reaction: Reaction) -> OrderedDict:
    """Convert a cobra Reaction object to a dictionary.

    Parameters
    ----------
    reaction : cobra.Reaction
        The cobra.Reaction to convert to dictionary.

    Returns
    -------
    dict
        The converted dictionary object.

    See Also
    --------
    _reaction_from_dict : Convert a dictionary to a cobra Reaction object.

    """
    new_reaction = OrderedDict()
    for key in _REQUIRED_REACTION_ATTRIBUTES:
        if key == "id":
            new_reaction[key] = _fix_type(F_REPLACE[F_REACTION_REV](reaction.id))
        else:
            new_reaction[key] = _fix_type(getattr(reaction, key))
    if F_REPLACE and F_SPECIE_REV in F_REPLACE:
        mets = {
            F_REPLACE[F_SPECIE_REV](met.id): stoic
            for met, stoic in reaction.metabolites.items()
        }
    else:
        mets = {met.id: stoic for met, stoic in reaction.metabolites.items()}
    new_reaction["metabolites"] = OrderedDict(mets)
    _update_optional(
        reaction,
        new_reaction,
        _OPTIONAL_REACTION_ATTRIBUTES,
        _ORDERED_OPTIONAL_REACTION_KEYS,
    )
    return new_reaction


def _reaction_from_dict(
    reaction: Dict, model: Model, f_replace: Dict = F_REPLACE  # noqa:    W0102
) -> Reaction:
    """Convert a dictionary to a cobra Reaction object.

    Parameters
    ----------
    reaction : dict
        The dictionary to convert to cobra.Reaction .
    model : cobra.Model
        The model to which the reaction should associate with.

    Returns
    -------
    cobra.Reaction
        The converted cobra.Reaction object.

    See Also
    --------
    _reaction_to_dict : Convert a cobra Reaction object to a dictionary.

    """
    if f_replace is None:
        f_replace = {}

    new_reaction = Reaction(_fix_id_from_dict(reaction["id"], "Reaction", f_replace))
    [
        setattr(new_reaction, k, _fix_value_from_dict(k, v))
        for k, v in reaction.items()
        if k
        not in {
            "id",
            "objective_coefficient",
            "reversibility",
            "reaction",
            "metabolites",
        }
    ]

    new_reaction.add_metabolites(
        OrderedDict(
            (
                model.metabolites.get_by_id(
                    _fix_id_from_dict(str(met), "Metabolite", f_replace)
                ),
                coeff,
            )
            for met, coeff in reaction.get("metabolites", {}).items()
        )
    )
    return new_reaction


def group_to_dict(group: "Group") -> Dict:
    new_group = OrderedDict()
    for key in _REQUIRED_GROUP_ATTRIBUTES:
        if key != "members":
            if key == "id":
                new_group[key] = _fix_type(F_REPLACE[F_GROUP_REV](group.id))
            else:
                new_group[key] = _fix_type(getattr(group, key))
                continue
        members = []
        for member in group.members:
            idRef = member.id
            if isinstance(member, Reaction):
                idRef = F_REPLACE[F_REACTION_REV](member.id)
            elif isinstance(member, Gene):
                idRef = F_REPLACE[F_GENE_REV](member.id)
            elif isinstance(member, Metabolite):
                idRef = F_REPLACE[F_SPECIE_REV](member.id)
            elif isinstance(member, Group):
                idRef = F_REPLACE[F_GROUP_REV](member.id)
            json_member = {"idRef": idRef, "type": type(member).__name__}
            members.append(json_member)
        new_group["members"] = members
    _update_optional(
        group, new_group, _OPTIONAL_GROUP_ATTRIBUTES, _ORDERED_OPTIONAL_GROUP_KEYS
    )
    return new_group


def group_from_dict(
    group: Dict, model: Model, f_replace=F_REPLACE # noqa:    W0102
) -> Group:
    if f_replace is None:
        f_replace = {}

    new_group = Group(_fix_id_from_dict(group["id"], "Group"))
    [
        setattr(new_group, k, _fix_value_from_dict(k, v))
        for k, v in group.items()
        if k not in {"id", "members"}
    ]
    cobra_members = [
        _get_by_id(
            _fix_id_from_dict(member["idRef"], member["type"], f_replace),
            member["type"],
            model,
        )
        for member in group["members"]
    ]
    new_group.add_members(cobra_members)
    return new_group


def model_to_dict(
    model: Model, sort: bool = False, f_replace: dict = F_REPLACE  # noqa:    W0102
) -> OrderedDict:
    """Convert a cobra Model to a dictionary.

    Parameters
    ----------
    model : cobra.Model
        The model to reformulate as a dict.
    sort : bool, optional
        Whether to sort the metabolites, reactions, and genes or maintain the
        order defined in the model (default False).
    f_replace : dict of replacement functions for id replacement
        Dictionary of replacement functions for gene, specie, and reaction.
        By default, the following id changes are performed on import:
        clip G_ from genes, clip M_ from species, clip R_ from reactions
        If no replacements should be performed, set f_replace={} or None

    Returns
    -------
    OrderedDict
        A dictionary with keys: 'genes', 'compartments', 'id',
        'metabolites', 'notes' and 'reactions'; where 'metabolites', 'genes'
        and 'metabolites' are in turn lists of dictionaries holding all
        attributes to form the corresponding object.

    See Also
    --------
    model_from_dict : Convert a dictionary to a cobra Model.

    """
    if f_replace is None:
        f_replace = {}

    obj = OrderedDict()
    obj["metabolites"] = [
        _metabolite_to_dict(met, f_replace) for met in model.metabolites
    ]
    obj["reactions"] = list(map(_reaction_to_dict, model.reactions))
    obj["genes"] = list(map(_gene_to_dict, model.genes))
    obj["groups"] = list(map(group_to_dict, model.groups))

    # sbml meta info
    if hasattr(model, "_sbml"):
        obj["sbml_info"] = OrderedDict({key: _fix_type(value)
                                        for key, value in model._sbml.items()})
    obj["id"] = model.id
    _update_optional(
        model, obj, _OPTIONAL_MODEL_ATTRIBUTES, _ORDERED_OPTIONAL_MODEL_KEYS
    )
    if sort:
        obj["metabolites"].sort(key=lambda x: x["id"])
        obj["reactions"].sort(key=lambda x: x["id"])
        obj["genes"].sort(key=lambda x: x["id"])
        obj["groups"].sort(key=lambda x: x["id"])
    return obj


def model_from_dict(obj: Dict) -> Model:
    """Build a cobra Model from a dictionary.

    Models stored in JSON are first formulated as a dictionary that can be read
    to a cobra Model using this function.

    Parameters
    ----------
    obj : dict
        A dictionary with keys: 'genes', 'compartments', 'id',
        'metabolites', 'notes' and 'reactions'; where 'metabolites', 'genes'
        and 'metabolites' are in turn lists of dictionaries holding all
        attributes to form the corresponding object.

    Returns
    -------
    cobra.Model
        The generated model.

    Raises
    ------
    ValueError
        If `obj` has no 'reactions' attribute.

    See Also
    --------
    model_to_dict : Convert a cobra Model to a dictionary.

    """
    if "reactions" not in obj:
        raise ValueError("Object has no .reactions attribute. Cannot load.")
    model = Model()
    model.add_metabolites(
        [_metabolite_from_dict(metabolite) for metabolite in obj["metabolites"]]
    )
    model.genes.extend([gene_from_dict(gene) for gene in obj["genes"]])
    model.add_reactions(
        [_reaction_from_dict(reaction, model) for reaction in obj["reactions"]]
    )
    objective_reactions = [
        rxn for rxn in obj["reactions"] if rxn.get("objective_coefficient", 0) != 0
    ]
    coefficients = {
        model.reactions.get_by_id(F_REPLACE[F_REACTION](rxn["id"])): rxn[
            "objective_coefficient"
        ]
        for rxn in objective_reactions
    }
    if "groups" in obj:
        model.add_groups([group_from_dict(group, model) for group in obj["groups"]])
    set_objective(model, coefficients)

    # sbml meta info
    if "sbml_info" in obj:
        meta = {k: _fix_value_from_dict(k, v) for k, v in obj["sbml_info"].items()}
        model._sbml = meta

    [
        setattr(model, k, _fix_value_from_dict(k, v))
        for k, v in obj.items()
        if k in {"id", "name", "compartments", "annotation", "notes"}
    ]

    return model
