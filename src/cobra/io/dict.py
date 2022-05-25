"""Provide functions for cobrapy objects to generic Python objects and vice-versa."""
import itertools
import re
from collections import OrderedDict, defaultdict
from functools import partial
from typing import (
    TYPE_CHECKING,
    Callable,
    Dict,
    List,
    Set,
    Tuple,
    Union,
)

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

_REACTION_DICT = {
    "id": "",
    "name": None,
    "lower_bound": None,
    "upper_bound": None,
    "subsystem": "",
    "notes": {},
    "annotation": {},
}

_METABOLITE_DICT = {
    "id": "",
    "name": None,
    "compartment": None,
    "charge": None,
    "formula": None,
    "_bound": 0,
    "notes": {},
    "annotation": {},
}

_GENE_DICT = {
    "id": "",
    "name": None,
    "notes": {},
    "annotation": {},
}

_GROUP_DICT = {
    "id": "",
    "name": "",
    "kind": "",
    "notes": {},
    "annotation": {},
}

_MODEL_DICT = {
    "id": "",
    "name": None,
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


def _f_replace_object_id(cobra_object, f_replace=F_REPLACE) -> str:  # noqa:    W0102
    if f_replace is None:
        f_replace = {}
    if f_replace:
        if isinstance(cobra_object, Reaction):
            return F_REPLACE[F_REACTION_REV](cobra_object.id)
        elif isinstance(cobra_object, Gene):
            return F_REPLACE[F_GENE_REV](cobra_object.id)
        elif isinstance(cobra_object, Metabolite):
            return F_REPLACE[F_SPECIE_REV](cobra_object.id)
        elif isinstance(cobra_object, Group):
            return F_REPLACE[F_GROUP_REV](cobra_object.id)
    return cobra_object.id


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


def _object_to_dict(
    cobra_object: "Object", optional_attribute_dict: Dict, _f_replace_function: Callable
) -> OrderedDict:
    """Update `new_dict` with optional attributes from `cobra_object`.

    Parameters
    ----------
    cobra_object : cobra.Object
        The cobra Object to update optional attributes from.
    optional_attribute_dict : dict
        The dictionary to use as default value lookup store.

    Raises
    ------
    IndexError
        If key in `ordered_keys` is not found in `optional_attribute_dict`.
    AttributeError
        If key in `ordered_keys` is not found in `cobra_object`.

    """
    optional_attribute_dict = optional_attribute_dict.copy()
    state = cobra_object.__getstate__()
    state["id"] = _f_replace_function(cobra_object)
    state_fixed = {
        (re.sub("^_", "", key) if key != "_bound" else key): state[key]
        for key in state.keys()
    }
    if (
        state_fixed["notes"] is None
        or state_fixed["notes"].notes_xhtml is None
        or len(state_fixed["notes"].notes_xhtml) == 0
    ):
        optional_attribute_dict.pop("notes")
    cobra_dict = OrderedDict(
        {
            key: _fix_type(state_fixed[key])
            for key, default in optional_attribute_dict.items()
            if state_fixed[key] is not None and state_fixed[key] != default
        }
    )
    return cobra_dict



def object_from_dict(new_object: "Object", object_dict: Dict):
    [
        setattr(new_object, k, _fix_value_from_dict(k, v))
        for k, v in object_dict.items()
        if k
        not in {
            "id",
            "objective_coefficient",
            "reversibility",
            "reaction",
            "metabolites",
            "members",
        }
    ]


def _metabolite_to_dict(
    metabolite: Metabolite, f_replace_func: Callable
) -> OrderedDict:
    """Convert a cobra Metabolite object to dictionary.

    Parameters
    ----------
    metabolite : cobra.Metabolite
        The cobra.Metabolite to convert to dictionary.
    f_replace_func : partial replacement function for id replacement
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
    return _object_to_dict(
        metabolite, _METABOLITE_DICT, _f_replace_function=f_replace_func
    )


def _metabolite_from_dict(
    metabolite: Dict, f_replace: dict = F_REPLACE  # noqa:    W0102
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
    object_from_dict(new_metabolite, metabolite)
    return new_metabolite


def _gene_to_dict(gene: Gene, f_replace_func: Callable) -> OrderedDict:
    """Convert a cobra Gene object to dictionary.

    Parameters
    ----------
    gene : cobra.Gene
        The cobra.Gene to convert to dictionary.
    f_replace_func : Callable
        partial replacement function for id replacement
        By default, the following id changes are performed on import:
        clip G_ from genes, clip M_ from species, clip R_ from reactions
        If no replacements should be performed, set f_replace={} or None

    Returns
    -------
    dict
        The converted dictionary object.

    See Also
    --------
    _gene_from_dict : Convert a dictionary to cobra Gene object.

    """
    return _object_to_dict(gene, _GENE_DICT, _f_replace_function=f_replace_func)


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
    object_from_dict(new_gene, gene)
    return new_gene


def _reaction_to_dict(reaction: Reaction, f_replace_func: Callable) -> OrderedDict:
    """Convert a cobra Reaction object to a dictionary.

    Parameters
    ----------
    reaction : cobra.Reaction
        The cobra.Reaction to convert to dictionary.
    f_replace_func : Callable
        partial replacement function for id replacement
        By default, the following id changes are performed on import:
        clip G_ from genes, clip M_ from species, clip R_ from reactions
        If no replacements should be performed, set f_replace={} or None

    Returns
    -------
    dict
        The converted dictionary object.

    See Also
    --------
    _reaction_from_dict : Convert a dictionary to a cobra Reaction object.

    """
    new_reaction = _object_to_dict(
        reaction, _REACTION_DICT, _f_replace_function=f_replace_func
    )
    new_reaction["metabolites"] = OrderedDict(
        {f_replace_func(met): stoic for met, stoic in reaction.metabolites.items()}
    )
    new_reaction["gene_reaction_rule"] = _fix_type(reaction.gene_reaction_rule)
    new_reaction["objective_coefficient"] = _fix_type(reaction.objective_coefficient)
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
    object_from_dict(new_reaction, reaction)

    new_reaction.add_metabolites(
        {
            model.metabolites.get_by_id(
                _fix_id_from_dict(str(met), "Metabolite", f_replace)
            ): coeff
            for met, coeff in reaction.get("metabolites", {}).items()
        }
    )
    return new_reaction


def group_to_dict(group: "Group", f_replace_func: Callable) -> Dict:
    new_group = _object_to_dict(group, _GROUP_DICT, _f_replace_function=f_replace_func)
    new_group["members"] = [
        {"idRef": f_replace_func(member), "type": type(member).__name__}
        for member in group.members
    ]
    return new_group


def group_from_dict(
    group: Dict, model: Model, f_replace=F_REPLACE  # noqa:    W0102
) -> Group:
    if f_replace is None:
        f_replace = {}

    new_group = Group(_fix_id_from_dict(group["id"], "Group"))
    object_from_dict(new_group, group)
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

    _f_replace_function = partial(_f_replace_object_id, f_replace=f_replace)

    obj = _object_to_dict(model, _MODEL_DICT, _f_replace_function)
    obj["metabolites"] = [
        _metabolite_to_dict(met, _f_replace_function) for met in model.metabolites
    ]
    obj["reactions"] = [
        _reaction_to_dict(rxn, _f_replace_function) for rxn in model.reactions
    ]
    obj["genes"] = [_gene_to_dict(gene, _f_replace_function) for gene in model.genes]
    obj["groups"] = [
        group_to_dict(group, _f_replace_function) for group in model.groups
    ]

    # sbml meta info
    if hasattr(model, "_sbml"):
        obj["sbml_info"] = OrderedDict(
            {key: _fix_type(value) for key, value in model._sbml.items()}
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
        model._sbml = {
            k: _fix_value_from_dict(k, v) for k, v in obj["sbml_info"].items()
        }

    [
        setattr(model, k, _fix_value_from_dict(k, v))
        for k, v in obj.items()
        if k in {"id", "name", "compartments", "annotation", "notes"}
    ]

    return model
