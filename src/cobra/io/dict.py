"""Provide functions for cobrapy objects to generic Python objects and vice-versa."""

from collections import OrderedDict
from operator import attrgetter, itemgetter
from typing import TYPE_CHECKING, Dict, List, Sequence, Set, Union

import numpy as np

from ..core import Gene, Metabolite, Model, Reaction
from ..util.solver import set_objective


if TYPE_CHECKING:
    from cobra import Object

_REQUIRED_REACTION_ATTRIBUTES = [
    "id",
    "name",
    "metabolites",
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

_ORDERED_OPTIONAL_MODEL_KEYS = ["name", "compartments", "notes", "annotation"]
_OPTIONAL_MODEL_ATTRIBUTES = {
    "name": None,
    #  "description": None, should not actually be included
    "compartments": [],
    "notes": {},
    "annotation": {},
}


def _fix_type(
    value: Union[str, np.float, np.bool, Set, Dict]
) -> Union[str, float, bool, List, OrderedDict]:
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
    if isinstance(value, np.float):
        return float(value)
    if isinstance(value, np.bool):
        return bool(value)
    if isinstance(value, set):
        return list(value)
    if isinstance(value, dict):
        return OrderedDict((key, value[key]) for key in sorted(value))
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
        if value is None or value == default:
            continue
        new_dict[key] = _fix_type(value)


def _metabolite_to_dict(metabolite: Metabolite) -> OrderedDict:
    """Convert a cobra Metabolite object to dictionary.

    Parameters
    ----------
    metabolite : cobra.Metabolite
        The cobra.Metabolite to convert to dictionary.

    Returns
    -------
    dict
        The converted dictionary object.

    See Also
    --------
    _metabolite_from_dict : Convert a dictionary to cobra Metabolite object.

    """
    new_metabolite = OrderedDict()
    for key in _REQUIRED_METABOLITE_ATTRIBUTES:
        new_metabolite[key] = _fix_type(getattr(metabolite, key))
    _update_optional(
        metabolite,
        new_metabolite,
        _OPTIONAL_METABOLITE_ATTRIBUTES,
        _ORDERED_OPTIONAL_METABOLITE_KEYS,
    )
    return new_metabolite


def _metabolite_from_dict(metabolite: Dict) -> Metabolite:
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
    new_metabolite = Metabolite()
    for k, v in metabolite.items():
        setattr(new_metabolite, k, v)
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
        new_gene[key] = _fix_type(getattr(gene, key))
    _update_optional(
        gene, new_gene, _OPTIONAL_GENE_ATTRIBUTES, _ORDERED_OPTIONAL_GENE_KEYS
    )
    return new_gene


def gene_from_dict(gene: Dict) -> Gene:
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
    new_gene = Gene(gene["id"])
    for k, v in gene.items():
        setattr(new_gene, k, v)
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
        if key != "metabolites":
            if key == "lower_bound" and (
                np.isnan(reaction.lower_bound) or np.isinf(reaction.lower_bound)
            ):
                new_reaction[key] = str(_fix_type(getattr(reaction, key)))
            elif key == "upper_bound" and (
                np.isnan(reaction.upper_bound) or np.isinf(reaction.upper_bound)
            ):
                new_reaction[key] = str(_fix_type(getattr(reaction, key)))
            else:
                new_reaction[key] = _fix_type(getattr(reaction, key))
            continue
        mets = OrderedDict()
        for met in sorted(reaction.metabolites, key=attrgetter("id")):
            mets[str(met)] = reaction.metabolites[met]
        new_reaction["metabolites"] = mets
    _update_optional(
        reaction,
        new_reaction,
        _OPTIONAL_REACTION_ATTRIBUTES,
        _ORDERED_OPTIONAL_REACTION_KEYS,
    )
    return new_reaction


def _reaction_from_dict(reaction: Dict, model: Model) -> Reaction:
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
    new_reaction = Reaction()
    for k, v in reaction.items():
        if k in {"objective_coefficient", "reversibility", "reaction"}:
            continue
        elif k == "metabolites":
            new_reaction.add_metabolites(
                OrderedDict(
                    (model.metabolites.get_by_id(str(met)), coeff)
                    for met, coeff in v.items()
                )
            )
        else:
            if k == "lower_bound" or k == "upper_bound":
                setattr(new_reaction, k, float(v))
            else:
                setattr(new_reaction, k, v)
    return new_reaction


def model_to_dict(model: Model, sort: bool = False) -> OrderedDict:
    """Convert a cobra Model to a dictionary.

    Parameters
    ----------
    model : cobra.Model
        The model to reformulate as a dict.
    sort : bool, optional
        Whether to sort the metabolites, reactions, and genes or maintain the
        order defined in the model (default False).

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
    obj = OrderedDict()
    obj["metabolites"] = list(map(_metabolite_to_dict, model.metabolites))
    obj["reactions"] = list(map(_reaction_to_dict, model.reactions))
    obj["genes"] = list(map(_gene_to_dict, model.genes))
    obj["id"] = model.id
    _update_optional(
        model, obj, _OPTIONAL_MODEL_ATTRIBUTES, _ORDERED_OPTIONAL_MODEL_KEYS
    )
    if sort:
        get_id = itemgetter("id")
        obj["metabolites"].sort(key=get_id)
        obj["reactions"].sort(key=get_id)
        obj["genes"].sort(key=get_id)
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
        model.reactions.get_by_id(rxn["id"]): rxn["objective_coefficient"]
        for rxn in objective_reactions
    }
    set_objective(model, coefficients)
    for k, v in obj.items():
        if k in {"id", "name", "notes", "compartments", "annotation"}:
            setattr(model, k, v)
    return model
