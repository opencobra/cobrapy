
## Phenotype Phase Plane

# This example is available as an IPython [notebook](http://nbviewer.ipython.or
# g/github/opencobra/cobrapy/blob/master/documentation_builder/phenotype_phase_
# plane.ipynb).

# Load iJO1366 as a test model and import cobra

from time import time

import cobra
from cobra.test import ecoli_pickle, create_test_model

from cobra.flux_analysis.phenotype_phase_plane import \
    calculate_phenotype_phase_plane

model = create_test_model(ecoli_pickle)
model
# Output:
# <Model iJO1366 at 0x5b0abd0>

# We want to make a phenotype phase plane to evaluate uptakes of Glucose and
# Oxygen.
# 
# With [matplotlib](http://matplotlib.org) installed, this is as simple as

data = calculate_phenotype_phase_plane(model, "EX_glc_e", "EX_o2_e")
data.plot_matplotlib();

# If [brewer2mpl](https://pypi.python.org/pypi/brewer2mpl/) is installed, other
# color schemes can be used as well

data.plot_matplotlib("Pastel1")
data.plot_matplotlib("Dark2");

# The number of points which are plotted in each dimension can also be changed

calculate_phenotype_phase_plane(model, "EX_glc_e", "EX_o2_e",
                                reaction1_npoints=20,
                                reaction2_npoints=20).plot_matplotlib();

# The code can also use multiple processes to speed up calculations

start_time = time()
calculate_phenotype_phase_plane(model, "EX_glc_e", "EX_o2_e", n_processes=1)
print("took %.2f seconds with 1 process" % (time() - start_time))
start_time = time()
calculate_phenotype_phase_plane(model, "EX_glc_e", "EX_o2_e", n_processes=4)
print("took %.2f seconds with 4 process" % (time() - start_time))
# Prints:
# took 5.97 seconds with 1 process
# took 2.97 seconds with 4 process
