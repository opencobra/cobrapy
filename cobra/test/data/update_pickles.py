#!/usr/bin/env python
# This script regenerates pickles of cobra Models.  Should be
# performed after updating core classes to prevent subtle bugs.
try:
    from cPickle import load, dump
except:
    from pickle import load, dump

from os.path import isfile

from cobra import Model
from cobra.version import get_version
from cobra.io import read_sbml_model, read_legacy_sbml, write_sbml_model
from cobra.io import save_matlab_model, save_json_model
from cobra.test import create_test_model

model_names = ['salmonella', 'iJO1366',]

for model_name in model_names:
    # read in old pickle and model from sbml
    model_pickle = model_name + '.pickle'
    if model_name == "iJO1366":
        new_model = read_legacy_sbml(model_name + '.xml')
    else:
        new_model = read_sbml_model(model_name + '.xml')
    # update other attributes
    if isfile(model_name + ".genes"):
        with open(model_name + ".genes", "rb") as infile:
            gene_names = load(infile)
        for gene in new_model.genes:
            gene.name = gene_names[gene.id]
    if isfile(model_name + ".media"):
        with open(model_name + ".media", "rb") as infile:
            new_model.media_compositions = load(infile)
    new_model._cobra_version = get_version()
    # write out new pickle
    with open(model_pickle, 'wb') as outfile:
        dump(new_model, outfile, protocol=2)
    # write out other formats for iJO1366
    if model_name == "iJO1366":
        save_matlab_model(new_model, model_name + ".mat")
        save_json_model(new_model, model_name + ".json")
    if model_name == "salmonella":
        write_sbml_model(new_model, model_name + "_fbc.xml")
