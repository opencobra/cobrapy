
## Simulating Deletions

# This example is available as an IPython [notebook](http://nbviewer.ipython.or
# g/github/opencobra/cobrapy/blob/master/documentation_builder/deletions.ipynb)
# .

from time import time


from cobra.test import create_test_model, salmonella_pickle, ecoli_pickle
from cobra.flux_analysis import single_deletion
from cobra.flux_analysis import double_deletion


cobra_model = create_test_model(salmonella_pickle)


### Single Deletions

                Perform all single gene deletions on a model
                
start = time()  # start timer()
growth_rates, statuses = single_deletion(cobra_model)
print("All single gene deletions completed in %.2f sec" % (time() - start))
# Prints:
# All single gene deletions completed in 4.01 sec

# These can also be done for only a subset of genes

single_deletion(cobra_model, element_list=cobra_model.genes[:100]);


# Single deletions can also be run on reactions

start = time()  # start timer()
growth_rates, statuses = single_deletion(cobra_model, element_type="reaction")
print("All single reaction deletions completed in %.2f sec" % (time() - start))
# Prints:
# All single reaction deletions completed in 7.41 sec

### Double Deletions

# Double deletions run in a similar way

start = time()  # start timer()
double_deletion(cobra_model, element_list_1=cobra_model.genes[:100])
print("Double gene deletions for 100 genes completed in %.2f sec" % (time() - start))
# Prints:
# Double gene deletions for 100 genes completed in 4.94 sec

# By default, the double deletion function will automatically use
# multiprocessing, splitting the task over up to 4 cores if they are available.
# The number of cores can be manually sepcified as well. Setting use of a
# single core will disable use of the multiprocessing library, which often aids
# debuggging.

start = time()  # start timer()
double_deletion(cobra_model, element_list_1=cobra_model.genes[:100],
                number_of_processes=2)
t1 = time() - start
print("Double gene deletions for 100 genes completed in %.2f sec with 2 cores" % t1)

start = time()  # start timer()
double_deletion(cobra_model, element_list_1=cobra_model.genes[:100],
                number_of_processes=1)
t2 = time() - start
print("Double gene deletions for 100 genes completed in %.2f sec with 1 core" % t2)

print("Speedup of %.2fx" % (t2/t1))
# Prints:
# Double gene deletions for 100 genes completed in 4.02 sec with 2 cores
# Double gene deletions for 100 genes completed in 6.77 sec with 1 core
# Speedup of 1.69x

# Double deletions can also be run for reactions

start = time()
double_deletion(cobra_model, element_list_1=cobra_model.reactions[:100],
                element_type="reaction")
t = time() - start
print("Double reaction deletions for 100 reactions completed in %.2f sec" % t)
# Prints:
# Double reaction deletions for 100 reactions completed in 0.93 sec

# If pandas is installed, the results can be returned formatted as a
# pandas.DataFrame

frame = double_deletion(cobra_model, element_list_1=cobra_model.reactions[300:308],
                        element_type="reaction", return_frame=True)
frame[frame < 1e-9] = 0.  # round small values to 0
frame
