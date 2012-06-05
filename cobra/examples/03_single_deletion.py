#cobra/examples/03_single_deletion.py
#
#This file provides a simple example of how to perform
#a single deletion simulation
from cobra.flux_analysis import single_deletion
from cPickle import load, dump
from time import time
from math import floor
from cobra.manipulation import initialize_growth_medium
from cobra.test import salmonella_pickle #This is the name of the test file
with open(salmonella_pickle) as in_file:
    cobra_model = load(in_file)


initialize_growth_medium(cobra_model, 'LB')
#Expected growth rates for the salmonella model with deletions in LB medium
the_loci =  ['STM4081', 'STM0247', 'STM3867', 'STM2952']
the_genes = tpiA, metN, atpA, eno = map(cobra_model.genes.get_by_id, the_loci)
id_to_name = dict([(x.id, x.name) for x in the_genes])
the_growth_rates = {tpiA.id:2.41, metN.id:2.43, atpA.id:1.87, eno.id:1.81} #expected growth rates after deletion

#Perform a single gene deletion
the_results = single_deletion(cobra_model, [tpiA])


gene_list = the_growth_rates.keys()
#Perform deletions for all genes in the list
start_time = time()

rates, statuses, problems = single_deletion(cobra_model,
                                            gene_list)
for the_gene, v in statuses.items():
    if v != 'optimal':
        print '\t\tdeletion %s was not optimal'%the_gene
for the_gene, v in rates.items():
    v = floor(100*v)/100
    if v != the_growth_rates[the_gene]:
        print '\t\tFAILED: %s simulation (%1.3f) != expectation (%1.2f)'%(id_to_name[the_gene],
                                                                              v,
                                                                              the_growth_rates[the_gene])
    else:
        print '\t\tPASSED: %s simulation (%1.3f) ~= expectation (%1.2f)'%(id_to_name[the_gene],
                                                                              v,
                                                                              the_growth_rates[the_gene])


print '\t\tsingle deletion time: %f seconds'%(time() - start_time)


