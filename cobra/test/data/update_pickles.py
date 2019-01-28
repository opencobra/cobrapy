#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import absolute_import

from collections import OrderedDict
from json import dump as json_dump

import cobra
from cobra.io import (
    load_matlab_model, read_sbml_model, save_json_model, save_matlab_model,
    write_sbml_model)
from cobra.io.sbml3 import write_sbml2


# This script regenerates pickles of cobra Models.  Should be
# performed after updating core classes to prevent subtle bugs.
try:
    from cPickle import load, dump
except:
    from pickle import load, dump


# ecoli
ecoli_model = read_sbml_model("iJO1366.xml")
with open("iJO1366.pickle", "wb") as outfile:
    dump(ecoli_model, outfile, protocol=2)

# salmonella
salmonella = read_sbml_model("salmonella.xml")
with open("salmonella.genes", "rb") as infile:
    gene_names = load(infile)
for gene in salmonella.genes:
    gene.name = gene_names[gene.id]
with open("salmonella.media", "rb") as infile:
    salmonella.media_compositions = load(infile)
with open("salmonella.pickle", "wb") as outfile:
    dump(salmonella, outfile, protocol=2)

# create mini model from textbook
textbook = read_sbml_model("textbook.xml.gz")
mini = cobra.Model("mini_textbook")
mini.compartments = textbook.compartments


for r in textbook.reactions:
    if r.id in ("GLCpts", "PGI", "PFK", "FBA", "TPI", "GAPD", "PGK", "PGM",
                "ENO", "PYK", "EX_glc__D_e", "EX_h_e", "H2Ot", "ATPM",
                "PIt2r"):
        mini.add_reaction(r.copy())
mini.reactions.ATPM.upper_bound = mini.reactions.PGI.upper_bound
mini.objective = ["PFK", "ATPM"]  # No biomass, 2 reactions

# add in some information from iJO1366
mini.add_reaction(ecoli_model.reactions.LDH_D.copy())
mini.add_reaction(ecoli_model.reactions.EX_lac__D_e.copy())
r = cobra.Reaction("D_LACt2")
mini.add_reaction(r)
r.gene_reaction_rule = ecoli_model.reactions.D__LACt2pp.gene_reaction_rule
r.reaction = ecoli_model.reactions.D__LACt2pp.reaction.replace("_p", "_e")
mini.reactions.GLCpts.gene_reaction_rule = \
    ecoli_model.reactions.GLCptspp.gene_reaction_rule

# adjust bounds
for i in ["ATPM", "D_LACt2", "EX_lac__D_e", "LDH_D"]:
    mini.reactions.get_by_id(i).upper_bound = mini.reactions.PGI.upper_bound
for i in ["D_LACt2", "LDH_D"]:
    mini.reactions.get_by_id(i).lower_bound = mini.reactions.PGI.lower_bound
# set names and annotation
for g in mini.genes:
    try:
        tg = textbook.genes.get_by_id(g.id)
    except KeyError:
        continue
    g.name = tg.name
    g.annotation = tg.annotation
mini.reactions.sort()
mini.genes.sort()
mini.metabolites.sort()
# output to various formats
with open("mini.pickle", "wb") as outfile:
    dump(mini, outfile, protocol=2)
save_matlab_model(mini, "mini.mat")
save_json_model(mini, "mini.json", pretty=True)
write_sbml_model(mini, "mini_fbc2.xml")
write_sbml_model(mini, "mini_fbc2.xml.bz2")
write_sbml_model(mini, "mini_fbc2.xml.gz")
write_sbml2(mini, "mini_fbc1.xml", use_fbc_package=True)
write_sbml_model(mini, "mini_cobra.xml", use_fbc_package=False)
raven = load_matlab_model("raven.mat")
with open("raven.pickle", "wb") as outfile:
    dump(raven, outfile, protocol=2)

# TODO:these need a reference solutions rather than circular solution checking!

# fva results
fva_result = cobra.flux_analysis.flux_variability_analysis(textbook)
clean_result = OrderedDict()
for key in sorted(fva_result):
    clean_result[key] = {k: round(v, 5) for k, v in fva_result[key].items()}
with open("textbook_fva.json", "w") as outfile:
    json_dump(clean_result, outfile)

# fva with pfba constraint
fva_result = cobra.flux_analysis.flux_variability_analysis(textbook,
                                                           pfba_factor=1.1)
clean_result = OrderedDict()
for key in sorted(fva_result):
    clean_result[key] = {k: round(v, 5) for k, v in fva_result[key].items()}
with open("textbook_pfba_fva.json", "w") as outfile:
    json_dump(clean_result, outfile)

# textbook solution
solution = cobra.flux_analysis.parsimonious.pfba(textbook)
with open('textbook_solution.pickle', 'wb') as f:
    dump(solution, f, protocol=2)
