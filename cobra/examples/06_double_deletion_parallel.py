#cobra/examples/06_double_deletion_parallel.py
#
#Advanced user example illustrating how to perform double deletion studies on a multicore system
#
from cobra.flux_analysis.essentiality import deletion_analysis
from cPickle import load, dump
from time import time
from copy import deepcopy
from sys import argv
from os import listdir, getcwd, devnull
from os.path import abspath
from cobra.manipulation import initialize_growth_medium


print '%s is still under development.'%argv[0]
if len(argv) < 4:
    print 'Script must be called with the number of processes, a cobra.Model pickle,' +\
          '\nand the output file name.  An optional parameter is the number of genes. e.g., ' +\
          'python 03_double_deletion_parallel.py 12 ../test/data/salmonella.pickle' +\
          ' salmonella_double.pickle 100 \n will split a double deletion study of' + \
          '\nthe first 100 genes in the genes attribute of the cobra.Model in'+\
          '\n in salmonella.pickle.  If you want to run for all genes then do not' +\
          '\nenter the number of genes'
    print repr(argv)
    exit()

number_of_processes = int(argv[1])
in_filename = argv[2]
out_filename = argv[3]
with open(in_filename) as in_file:
    cobra_model = load(in_file)
if hasattr(cobra_model, 'media_compositions'):
    #If your pickle has a growth media dictionary then you can set the
    #medium component sources here.
    
    initialize_growth_medium(cobra_model, 'LB')

#When specifying the genes to delete use the locus ids because cobra.Gene objects
#are not thread safe and copy times might be excessive.
if len(argv) == 5:
    number_of_genes = int(argv[4])
    gene_list = [x.id for x in cobra_model.genes[:number_of_genes]]
else:
    gene_list = [x.id for x in cobra_model.genes]


start_time = time()
print 'running double deletion for %i genes on %i cores'%(len(gene_list), number_of_processes)
the_results = deletion_analysis(cobra_model, element_list=gene_list,n_processes=number_of_processes, deletion_type='double')
print 'took %1.2f seconds to do double deletion for %i genes'%(time() - start_time, len(gene_list))
with open(out_filename, 'w') as out_file:
    dump(the_results, out_file)
    print 'saved the results to %s'%out_filename


