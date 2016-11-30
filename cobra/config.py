from __future__ import absolute_import, print_function

import logging

from .parallel import SequentialView

logging.getLogger().setLevel(logging.ERROR)

log = logging.getLogger(__name__)

non_zero_flux_threshold = 1e-6
ndecimals = 6

# Determine available solver interfaces
solvers = {}

try:
    from optlang import glpk_interface

    solvers['glpk'] = glpk_interface
except ImportError:
    pass
try:
    from optlang import cplex_interface

    solvers['cplex'] = cplex_interface
except ImportError:
    pass

# Set default parallelization view
default_view = SequentialView()
