# This example demonstrates a single gene deletion simulation
from time import time

from cobra.flux_analysis.single_deletion import single_deletion
from cobra.manipulation import initialize_growth_medium
from cobra.test import create_test_model, salmonella_pickle  # test filename

cobra_model = create_test_model(salmonella_pickle)
initialize_growth_medium(cobra_model, 'LB')

target_genes = ['STM4081', 'STM0247', 'STM3867', 'STM2952']
# Expected growth rates for the salmonella model after a deletions in LB medium
expected_growth_rates = {
    "STM4081": 2.41,
    "STM0247": 2.43,
    "STM3867": 1.87,
    "STM2952": 1.81}

start_time = time()  # start timer

# Perform deletions for all genes in the list
rates, statuses = single_deletion(cobra_model, target_genes)

total_time = time() - start_time  # stop timer

# print out results
passed_string = 'PASSED: %s simulation (%1.3f) ~= expectation (%1.2f)'
failed_string = 'FAILED: %s simulation (%1.3f) != expectation (%1.2f)'
for gene_locus, rate in rates.items():
    # get gene name from gene locus (i.e. STM4081 -> tpiA)
    name = cobra_model.genes.get_by_id(gene_locus).name
    # test if the simulation failed
    if statuses[gene_locus] != "optimal":
        print("deletion failed for %s (%s)" % (name, gene_locus))
    if abs(rate - expected_growth_rates[gene_locus]) > 0.01:
        print(failed_string % (name, rate, expected_growth_rates[gene_locus]))
    else:
        print(passed_string % (name, rate, expected_growth_rates[gene_locus]))
print('single deletion time: %f seconds' % (total_time))
