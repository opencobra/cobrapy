#cobra/examples/04_change_objective_example.py
#
#This file shows how to change the targeted objective
#for the optimization function.
#
from cobra.flux_analysis import single_deletion
from cPickle import load, dump
from time import time
from math import floor

#Load in the example model file
test_directory = 'files/'
from cobra.manipulation import initialize_growth_medium
with open(test_directory + 'salmonella.pickle') as in_file:
    cobra_model = load(in_file)


cobra_model.optimize()
print 'solution for old objective (should be approximately 0.380):'
print cobra_model.solution.f

my_objective = cobra_model.reactions.get_by_id('3OAS140')
print 'Changing objective to %s'%my_objective.id
cobra_model.optimize(new_objective=my_objective)
print 'solution for %s (should be approximately 0.845):'%my_objective.id
print cobra_model.solution.f

