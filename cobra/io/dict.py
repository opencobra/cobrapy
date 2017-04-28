# -*- coding: utf-8 -*-

from __future__ import absolute_import

from collections import OrderedDict

from numpy import bool_, float_
from six import iteritems, string_types

from cobra.core import Gene, Metabolite, Model, Reaction
from cobra.util.solver import set_objective

_ORDERED_REACTION_KEYS = [
    "id", "name", "metabolites", "lower_bound", "upper_bound",
    "gene_reaction_rule"]
_REQUIRED_REACTION_ATTRIBUTES = {"id", "name", "metabolites", "lower_bound",
                                 "upper_bound", "gene_reaction_rule"}
_ORDERED_OPTIONAL_REACTION_KEYS = [
    "objective_coefficient", "variable_kind", "subsystem", "notes",
    "annotation"]
_OPTIONAL_REACTION_ATTRIBUTES = {
    "objective_coefficient": 0,
    "variable_kind": "continuous",
    "subsystem": "",
    "notes": {},
    "annotation": {},
}

_ORDERED_METABOLITE_KEYS = [
    "id", "name", "compartment"]
_REQUIRED_METABOLITE_ATTRIBUTES = {"id", "name", "compartment"}
_ORDERED_OPTIONAL_METABOLITE_KEYS = [
    "charge", "formula", "_bound", "_constraint_sense", "notes", "annotation"]
_OPTIONAL_METABOLITE_ATTRIBUTES = {
    "charge": None,
    "formula": None,
    "_bound": 0,
    "_constraint_sense": "E",
    "notes": {},
    "annotation": {},
}

_ORDERED_GENE_KEYS = ["id", "name"]
_REQUIRED_GENE_ATTRIBUTES = {"id", "name"}
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


def _fix_type(value):
    """convert possible types to str, float, and bool"""
    # Because numpy floats can not be pickled to json
    if isinstance(value, string_types):
        return str(value)
    if isinstance(value, float_):
        return float(value)
    if isinstance(value, bool_):
        return bool(value)
    if isinstance(value, set):
        return list(value)
    # handle legacy Formula type
    if value.__class__.__name__ == "Formula":
        return str(value)
    if value is None:
        return ''
    return value


def _update_optional(cobra_object, new_dict, optional_attribute_dict,
                     ordered_keys=None):
    """update new_dict with optional attributes from cobra_object"""
    if ordered_keys is not None:
        items = ((key, optional_attribute_dict[key]) for key in ordered_keys)
    else:
        items = iteritems(optional_attribute_dict)
    for key, default in items:
        value = getattr(cobra_object, key)
        if value is None or value == default:
            continue
        new_dict[key] = _fix_type(value)


def metabolite_to_dict(metabolite, ordered=False):
    if ordered:
        new_met = OrderedDict()
        for key in _ORDERED_METABOLITE_KEYS:
            new_met[key] = _fix_type(getattr(metabolite, key))
        _update_optional(metabolite, new_met, _OPTIONAL_METABOLITE_ATTRIBUTES,
                         _ORDERED_OPTIONAL_METABOLITE_KEYS)
    else:
        new_met = {key: _fix_type(getattr(metabolite, key))
                   for key in _REQUIRED_METABOLITE_ATTRIBUTES}
        _update_optional(metabolite, new_met, _OPTIONAL_METABOLITE_ATTRIBUTES)
    return new_met


def metabolite_from_dict(metabolite):
    new_metabolite = Metabolite()
    for k, v in iteritems(metabolite):
        setattr(new_metabolite, k, v)
    return new_metabolite


def gene_to_dict(gene, ordered=False):
    if ordered:
        new_gene = OrderedDict()
        for key in _ORDERED_GENE_KEYS:
            new_gene[key] = str(getattr(gene, key))
        _update_optional(gene, new_gene, _OPTIONAL_GENE_ATTRIBUTES,
                         _ORDERED_OPTIONAL_GENE_KEYS)
    else:
        new_gene = {key: str(getattr(gene, key))
                    for key in _REQUIRED_GENE_ATTRIBUTES}
        _update_optional(gene, new_gene, _OPTIONAL_GENE_ATTRIBUTES)
    return new_gene


def gene_from_dict(gene):
    new_gene = Gene(gene["id"])
    for k, v in iteritems(gene):
        setattr(new_gene, k, v)
    return new_gene


def reaction_to_dict(reaction, ordered=False):
    if ordered:
        new_reaction = OrderedDict()
        for key in _ORDERED_REACTION_KEYS:
            if key != "metabolites":
                new_reaction[key] = _fix_type(getattr(reaction, key))
                continue
            mets = OrderedDict()
            for met, coeff in iteritems(reaction.metabolites):
                mets[str(met)] = coeff
            new_reaction["metabolites"] = mets
        _update_optional(reaction, new_reaction, _OPTIONAL_REACTION_ATTRIBUTES,
                         _ORDERED_OPTIONAL_REACTION_KEYS)
    else:
        new_reaction = {key: _fix_type(getattr(reaction, key))
                        for key in _REQUIRED_REACTION_ATTRIBUTES
                        if key != "metabolites"}
        _update_optional(reaction, new_reaction, _OPTIONAL_REACTION_ATTRIBUTES)
        # set metabolites
        mets = {str(met): coeff for met, coeff
                in iteritems(reaction.metabolites)}
        new_reaction["metabolites"] = mets
    return new_reaction


def reaction_from_dict(reaction, model):
    new_reaction = Reaction()
    for k, v in iteritems(reaction):
        if k in {'objective_coefficient', 'reversibility', 'reaction'}:
            continue
        elif k == 'metabolites':
            new_reaction.add_metabolites(
                {model.metabolites.get_by_id(str(met)): coeff
                 for met, coeff in iteritems(v)})
        else:
            setattr(new_reaction, k, v)
    return new_reaction


def model_to_dict(model, ordered=False):
    """Convert model to a dict.

    Parameters
    ----------
    model : cobra.Model
        The model to reformulate as a dict
    ordered : bool, optional
        Whether to use ordered keys as determined by the JSON spec and
        elements being sorted by ID (default is arbitrary order).

    Returns
    -------
    dict
        A dictionary with elements, 'genes', 'compartments', 'id',
        'metabolites', 'notes' and 'reactions'; where 'metabolites', 'genes'
        and 'metabolites' are in turn lists with dictionaries holding all
        attributes to form the corresponding object.

    See Also
    --------
    cobra.io.model_from_dict
    """
    if ordered:
        # Currently rely on model.reactions and others having a specific order.
        # May need to sort by ID if this is not the case.
        obj = OrderedDict()
        obj["reactions"] = [
            reaction_to_dict(reaction, True) for reaction in model.reactions]
        obj["metabolites"] = [
            metabolite_to_dict(metabolite, True)
            for metabolite in model.metabolites]
        obj["genes"] = [gene_to_dict(gene, True) for gene in model.genes]
        obj["id"] = model.id
        _update_optional(model, obj, _OPTIONAL_MODEL_ATTRIBUTES,
                         _ORDERED_OPTIONAL_MODEL_KEYS)
    else:
        obj = dict(
            reactions=[reaction_to_dict(reaction)
                       for reaction in model.reactions],
            metabolites=[metabolite_to_dict(metabolite)
                         for metabolite in model.metabolites],
            genes=[gene_to_dict(gene) for gene in model.genes],
            id=model.id
        )
        _update_optional(model, obj, _OPTIONAL_MODEL_ATTRIBUTES)
    # add in the YAML version
    obj["version"] = 1
    return obj


def model_from_dict(obj):
    """Build a model from a dict.

    Models stored in json are first formulated as a dict that can be read to
    cobra model using this function.

    Parameters
    ----------
    obj : dict
        A dictionary with elements, 'genes', 'compartments', 'id',
        'metabolites', 'notes' and 'reactions'; where 'metabolites', 'genes'
        and 'metabolites' are in turn lists with dictionaries holding all
        attributes to form the corresponding object.

    Returns
    -------
    cora.core.Model
        The generated model.

    See Also
    --------
    cobra.io.model_to_dict
    """
    if 'reactions' not in obj:
        raise Exception('JSON object has no reactions attribute. Cannot load.')
    model = Model()
    model.add_metabolites(
        [metabolite_from_dict(metabolite) for metabolite in obj['metabolites']]
    )
    model.genes.extend([gene_from_dict(gene) for gene in obj['genes']])
    model.add_reactions(
        [reaction_from_dict(reaction, model) for reaction in obj['reactions']]
    )
    objective_reactions = [rxn for rxn in obj['reactions'] if
                           rxn.get('objective_coefficient', 0) != 0]
    coefficients = {
        model.reactions.get_by_id(rxn['id']): rxn['objective_coefficient'] for
        rxn in objective_reactions}
    set_objective(model, coefficients)
    for k, v in iteritems(obj):
        if k in {'id', 'name', 'notes', 'compartments', 'annotation'}:
            setattr(model, k, v)
    return model
