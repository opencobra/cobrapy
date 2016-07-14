from __future__ import absolute_import

import json
from warnings import warn

from .. import Model, Metabolite, Reaction, Gene
from six import iteritems, string_types

# Detect numpy types to replace them.
try:
    from numpy import float_, bool_
except ImportError:
    class float_:
        pass

    class bool_:
        pass

_REQUIRED_REACTION_ATTRIBUTES = {"id", "name", "metabolites", "lower_bound",
                                 "upper_bound", "gene_reaction_rule"}
_OPTIONAL_REACTION_ATTRIBUTES = {
    "objective_coefficient": 0,
    "variable_kind": "continuous",
    "subsystem": "",
    "notes": {},
    "annotation": {},
}

_REQUIRED_METABOLITE_ATTRIBUTES = {"id", "name", "compartment"}
_OPTIONAL_METABOLITE_ATTRIBUTES = {
    "charge": None,
    "formula": None,
    "_bound": 0,
    "_constraint_sense": "E",
    "notes": {},
    "annotation": {},
}

_REQUIRED_GENE_ATTRIBUTES = {"id", "name"}
_OPTIONAL_GENE_ATTRIBUTES = {
    "notes": {},
    "annotation": {},
}

_OPTIONAL_MODEL_ATTRIBUTES = {
    "name": None,
    #  "description": None, should not actually be included
    "compartments": {},
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
    # handle legacy Formula type
    if value.__class__.__name__ == "Formula":
        return str(value)
    if value is None:
        return ''
    return value


def _from_dict(obj):
    """build a model from a dict"""
    if 'reactions' not in obj:
        raise Exception('JSON object has no reactions attribute. Cannot load.')
    model = Model()
    # add metabolites
    new_metabolites = []
    for metabolite in obj['metabolites']:
        new_metabolite = Metabolite()
        for k, v in iteritems(metabolite):
            setattr(new_metabolite, k, v)
        new_metabolites.append(new_metabolite)
    model.add_metabolites(new_metabolites)
    # add genes
    for gene in obj['genes']:
        new_gene = Gene(gene["id"])
        for k, v in iteritems(gene):
            setattr(new_gene, k, v)
        model.genes.append(new_gene)
    # add reactions
    new_reactions = []
    for reaction in obj['reactions']:
        new_reaction = Reaction()
        for k, v in iteritems(reaction):
            if k == 'reversibility' or k == "reaction":
                continue
            elif k == 'metabolites':
                new_reaction.add_metabolites(
                    {model.metabolites.get_by_id(str(met)): coeff
                     for met, coeff in iteritems(v)})
            else:
                setattr(new_reaction, k, v)
        new_reactions.append(new_reaction)
    model.add_reactions(new_reactions)
    for k, v in iteritems(obj):
        if k in {'id', 'name', 'notes', 'compartments', 'annotation'}:
            setattr(model, k, v)
    return model


def _update_optional(cobra_object, new_dict, optional_attribute_dict):
    """update new_dict with optional attributes from cobra_object"""
    for key, default_value in iteritems(optional_attribute_dict):
        value = getattr(cobra_object, key)
        if value is not None and value != default_value:
            new_dict[key] = _fix_type(value)


def _to_dict(model):
    """convert the model to a dict"""
    new_reactions = []
    new_metabolites = []
    new_genes = []
    for reaction in model.reactions:
        new_reaction = {key: _fix_type(getattr(reaction, key))
                        for key in _REQUIRED_REACTION_ATTRIBUTES
                        if key != "metabolites"}
        _update_optional(reaction, new_reaction, _OPTIONAL_REACTION_ATTRIBUTES)
        # set metabolites
        mets = {str(met): coeff for met, coeff
                in iteritems(reaction._metabolites)}
        new_reaction['metabolites'] = mets
        new_reactions.append(new_reaction)
    for metabolite in model.metabolites:
        new_metabolite = {key: _fix_type(getattr(metabolite, key))
                          for key in _REQUIRED_METABOLITE_ATTRIBUTES}
        _update_optional(metabolite, new_metabolite,
                         _OPTIONAL_METABOLITE_ATTRIBUTES)
        new_metabolites.append(new_metabolite)
    for gene in model.genes:
        new_gene = {key: str(getattr(gene, key))
                    for key in _REQUIRED_GENE_ATTRIBUTES}
        _update_optional(gene, new_gene, _OPTIONAL_GENE_ATTRIBUTES)
        new_genes.append(new_gene)
    obj = {'reactions': new_reactions,
           'metabolites': new_metabolites,
           'genes': new_genes,
           'id': model.id,
           }

    _update_optional(model, obj, _OPTIONAL_MODEL_ATTRIBUTES)
    # add in the JSON version
    obj["version"] = 1
    return obj


def to_json(model):
    """Save the cobra model as a json string"""
    return json.dumps(_to_dict(model), allow_nan=False)


def from_json(jsons):
    """Load cobra model from a json string"""
    return _from_dict(json.loads(jsons))


def load_json_model(file_name):
    """Load a cobra model stored as a json file

    file_name : str or file-like object

    """
    # open the file
    should_close = False
    if isinstance(file_name, string_types):
        file_name = open(file_name, 'r')
        should_close = True

    model = _from_dict(json.load(file_name))

    if should_close:
        file_name.close()

    return model


def save_json_model(model, file_name, pretty=False):
    """Save the cobra model as a json file.

    model : :class:`~cobra.core.Model.Model` object

    file_name : str or file-like object

    """
    # open the file
    should_close = False
    if isinstance(file_name, string_types):
        file_name = open(file_name, 'w')
        should_close = True

    if pretty:
        dump_opts = {"indent": 4, "separators": (",", ": "), "sort_keys": True}
    else:
        dump_opts = {}

    json.dump(_to_dict(model), file_name, allow_nan=False, **dump_opts)

    if should_close:
        file_name.close()


json_schema = {
    "$schema": "http://json-schema.org/draft-04/schema#",
    "title": "COBRA",
    "description": "JSON representation of COBRA model",
    "type": "object",
    "properties": {
        "id": {"type": "string"},
        "name": {"type": "string"},
        "description": {"type": "string"},
        "version": {
            "type": "integer",
            "default": 1,
        },

        "reactions": {
            "type": "array",
            "items": {
                "type": "object",
                "properties": {
                    "id": {"type": "string"},
                    "name": {"type": "string"},
                    "metabolites": {
                        "type": "object",
                        "patternProperties": {
                            ".*": {"type": "number"},
                        }
                    },
                    "gene_reaction_rule": {"type": "string"},
                    "lower_bound": {"type": "number"},
                    "upper_bound": {"type": "number"},
                    "objective_coefficient": {
                        "type": "number",
                        "default": 0,
                    },
                    "variable_kind": {
                        "type": "string",
                        "pattern": "integer|continuous",
                        "default": "continuous"
                    },
                    "subsystem": {"type": "string"},
                    "notes": {"type": "object"},
                    "annotation": {"type": "object"},
                },
                "required": ["id", "name", "metabolites", "lower_bound",
                             "upper_bound", "gene_reaction_rule"],
                "additionalProperties": False,
            }
        },
        "metabolites": {
            "type": "array",
            "items": {
                "type": "object",
                "properties": {
                    "id": {"type": "string"},
                    "name": {"type": "string"},
                    "compartment": {
                        "type": "string",
                        "pattern": "[a-z]{1,2}"
                    },
                    "charge": {"type": "integer"},
                    "formula": {"type": "string"},
                    "_bound": {
                        "type": "number",
                        "default": 0
                    },
                    "_constraint_sense": {
                        "type": "string",
                        "default": "E",
                        "pattern": "E|L|G",
                    },
                    "notes": {"type": "object"},
                    "annotation": {"type": "object"},
                },
                "required": ["id", "name", "compartment"],
                "additionalProperties": False,
            }

        },
        "genes": {
            "type": "array",
            "items": {
                "type": "object",
                "properties": {
                    "id": {"type": "string"},
                    "name": {"type": "string"},
                    "notes": {"type": "object"},
                    "annotation": {"type": "object"},
                },
                "required": ["id", "name"],
                "additionalProperties": False,
            }

        },
        "compartments": {
            "type": "object",
            "patternProperties": {
                "[a-z]{1,2}": {"type": "string"}
            }
        },
        "notes": {"type": "object"},
        "annotation": {"type": "object"},
    },
    "required": ["id", "reactions", "metabolites", "genes"],
    "additionalProperties": False,
}
