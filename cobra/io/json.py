from __future__ import absolute_import

import json
from warnings import warn

from .. import Model, Metabolite, Reaction, Formula

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
    for metabolite_id, metabolite in obj['metabolites'].iteritems():
        new_metabolite = Metabolite()
        new_metabolite.id = metabolite_id
        for k, v in metabolite.iteritems():
            setattr(new_metabolite, k, v)
            # TODO test that these are the right attributes
        new_metabolites.append(new_metabolite)
    model.add_metabolites(new_metabolites)
    # add reactions
    new_reactions = []
    for reaction_id, reaction in obj['reactions'].iteritems():
        new_reaction = Reaction()
        new_reaction.id = reaction_id
        for k, v in reaction.iteritems():
            if k=='reversibility': continue
            # TODO test that these are the right attributes
            if k=='metabolites':
                for met, coeff in v.iteritems():
                    new_reaction.add_metabolites({model.metabolites.get_by_id(met): coeff})
                continue
            setattr(new_reaction, k, v)
        new_reactions.append(new_reaction)
    model.add_reactions(new_reactions)
    for k, v in obj.iteritems():
        if k in ['id', 'description', 'notes']:
            setattr(model, k, v)
    return model

def _to_dict(model, exclude_attributes=[]):
    new_reactions = {}; new_metabolites = {}
    for reaction in model.reactions:
        new_reaction = {}
        # set metabolites
        mets = {str(met): coeff for met, coeff in reaction._metabolites.iteritems()}
        new_reaction['metabolites'] = mets
        for key in ['name', 'reversibility', 'subsystem', 'upper_bound', 'lower_bound',
                    'objective_coefficient', 'notes']:
            if key in exclude_attributes: continue
            new_reaction[key] = getattr(reaction, key)
        new_reactions[reaction.id] = new_reaction
    for metabolite in model.metabolites:
        new_metabolite = {}
        for key in ['annotation', 'charge', 'compartment', 'formula', 'name', 'notes']:
            if key in exclude_attributes: continue
            new_metabolite[key] = str(getattr(metabolite, key))
        new_metabolites[metabolite.id] = new_metabolite
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
    if isinstance(file_name, str):
        file_name = open(file_name, 'w')

    # iterator
    json.dump(_to_dict(model, exclude_attributes=exclude_attributes), file_name)
