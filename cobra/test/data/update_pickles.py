#!/usr/bin/env python
#This script regenerates pickles of cobra Models.  Should be
#performed after updating core classes to prevent subtle bugs.
from cPickle import load, dump
from cobra import Model
from cobra.version import get_version
from cobra.io import read_sbml_model, read_legacy_sbml
from cobra.test import create_test_model

model_names = ['salmonella', 'iJO1366', 'Yersinia_pestis_CO92_iPC815']

for model_name in model_names:
    # read in old pickle and model from sbml
    model_pickle = model_name + '.pickle'
    old_model = create_test_model(model_pickle)
    if model_name == "iJO1366":
        new_model = read_legacy_sbml(model_name + '.xml')
    else:
        new_model = read_sbml_model(model_name + '.xml')
    # update any attributes from the old pickle
    for x in new_model.genes:    
        x.name = old_model.genes.get_by_id(x.id).name
    if hasattr(old_model, 'media_compositions'):
        new_model.media_compositions = old_model.media_compositions
    new_model._cobra_version = get_version()
    # write out new pickle
    with open(model_pickle, 'wb') as outfile:
        dump(new_model, outfile, protocol=2)
