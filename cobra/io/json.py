from __future__ import absolute_import

import json
from warnings import warn

from .. import Model, Metabolite, Reaction, Formula
from ..external.six import iteritems

def load_json_model(infile_path, variable_name=None):
    """Load a cobra model stored as a json file

    Parameters
    ----------
    infile_path : str

    """
    model = Model()
    with open(infile_path, 'r') as f:
        obj = json.load(f)
        
    if not 'reactions' in obj:
        raise Exception('JSON object has no reactions attribute. Cannot load.')
    # add metabolites
    new_metabolites = []
    for metabolite in obj['metabolites']:
        new_metabolite = Metabolite()
        for k, v in iteritems(metabolite):
            setattr(new_metabolite, k, v)
            # TODO test that these are the right attributes
        new_metabolites.append(new_metabolite)
    model.add_metabolites(new_metabolites)
    # add reactions
    new_reactions = []
    for reaction in obj['reactions']:
        new_reaction = Reaction()
        for k, v in iteritems(reaction):
            if k == 'reversibility' or k == "reaction":
                continue
            # TODO test that these are the right attributes
            elif k=='metabolites':
                new_reaction.add_metabolites(
                    {model.metabolites.get_by_id(met): coeff
                     for met, coeff in iteritems(v)})
            else:
                setattr(new_reaction, k, v)
        new_reactions.append(new_reaction)
    model.add_reactions(new_reactions)
    for k, v in iteritems(obj):
        if k in ['id', 'description', 'notes']:
            setattr(model, k, v)
    return model

_DEFAULT_REACTION_ATTRIBUTES = {
    'id', 'name', 'reversibility', 'subsystem', 'lower_bound', 'upper_bound',
    'objective_coefficient', 'notes', 'gene_reaction_rule', 'reaction'}

_DEFAULT_METABOLITE_ATTRIBUTES = {
    'id', 'annotation', 'charge', 'compartment', 'formula', 'name', 'notes'}

def _to_dict(model, exclude_attributes=[]):
    exclude_attributes = set(exclude_attributes)
    reaction_attributes = _DEFAULT_REACTION_ATTRIBUTES - exclude_attributes
    metabolite_attributes = _DEFAULT_METABOLITE_ATTRIBUTES - exclude_attributes
    new_reactions = []
    new_metabolites = []
    for reaction in model.reactions:
        new_reaction = {key: getattr(reaction, key)
                        for key in reaction_attributes}
        # set metabolites
        mets = {str(met): coeff for met, coeff in iteritems(reaction._metabolites)}
        new_reaction['metabolites'] = mets
        new_reactions.append(new_reaction)
    for metabolite in model.metabolites:
        new_metabolite = {key: str(getattr(metabolite, key))
                          for key in metabolite_attributes}
        new_metabolites.append(new_metabolite)
    obj = {'reactions': new_reactions,
           'metabolites': new_metabolites,
           'id': model.id,
           'description': model.description,
           'notes': model.notes}
    return obj

def to_json(model, exclude_attributes=[]):
    return json.dumps(_to_dict(model, exclude_attributes=exclude_attributes))

def save_json_model(model, file_name, exclude_attributes=[]):
    """Save the cobra model as a json file.

    Parameters
    ----------
    model : cobra.Model
    file_name : str or file-like object
    exclude_attributes : A list of reaction or metabolite attributes to ignore.
    Warning: ignoring attributes will make it impossible to reload the model with
    cobra.io.json.load_json_model()

    """
    # open the file
    should_close = False
    if isinstance(file_name, str):
        file_name = open(file_name, 'w')
        should_close = True

    json.dump(_to_dict(model, exclude_attributes=exclude_attributes), file_name)

    if should_close:
        file_name.close()
