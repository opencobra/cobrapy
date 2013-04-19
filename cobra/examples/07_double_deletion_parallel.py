#cobra/examples/07_double_deletion_parallel.py
#
#Advanced user example illustrating how to perform double deletion studies on a multicore system
#
from cobra.flux_analysis.double_deletion import double_deletion
from cPickle import load, dump
from time import time
from cobra.test import create_test_model

number_of_processes = 4 #Number of parallel processes to start
number_of_genes = 100 #Total number of genes to perform double deletion on
out_filename = 'double_deletion_results.pickle'

cobra_model = create_test_model()

#When specifying the genes to delete use the locus ids because cobra.Gene objects
#are not thread safe and copy times might be excessive.
gene_list = [x.id for x in cobra_model.genes[:number_of_genes]]


start_time = time()
print 'running double deletion for %i genes on %i cores'%(len(gene_list), number_of_processes)
the_results = double_deletion(cobra_model, element_list_1=gene_list,number_of_processes=number_of_processes)
print 'took %1.2f seconds to do double deletion for %i genes'%(time() - start_time, len(gene_list))
with open(out_filename, 'w') as out_file:
    dump(the_results, out_file)
    print 'saved the results to %s'%out_filename


