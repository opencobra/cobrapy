from __future__ import absolute_import, print_function

import logging

logging.getLogger().setLevel(logging.ERROR)

log = logging.getLogger(__name__)

non_zero_flux_threshold = 1e-6
ndecimals = 6
