#This script regenerates pickles of cobra Models.  Should be
#performed after updating core classes to prevent subtle bugs.
from cPickle import load, dump
from cobra import Model
from cobra.version import get_version
from cobra.io import read_sbml_model
from cobra.test import create_test_model
model_names = ['salmonella', 'iJO1366', 'Yersinia_pestis_CO92_iPC815']
for model_name in model_names:
    model_pickle = model_name + '.pickle'
    old_model = create_test_model(model_pickle)
    new_model = Model(old_model)
    new_model._cobra_version = get_version()
    dump(new_model, open(model_pickle, 'w'))
