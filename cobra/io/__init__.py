# -*- coding: utf-8 -*-

from __future__ import absolute_import

from warnings import warn

from cobra.io.json import load_json_model, save_json_model, to_json
from cobra.io.sbml3 import read_sbml_model, write_sbml_model

# These functions have other dependencies
try:
    import libsbml
    from cobra.io.sbml import read_legacy_sbml
    from cobra.io.sbml import write_cobra_model_to_sbml_file as\
        write_legacy_sbml
except ImportError:
    warn("cobra.io.sbml requires libsbml")
    libsbml = None
    read_legacy_sbml = None
    write_legacy_sbml = None

try:
    import scipy
    from cobra.io.mat import load_matlab_model, save_matlab_model
except ImportError:
    warn("cobra.io.mat requires scipy")
    scipy = None
    load_matlab_model = None
    save_matlab_model = None

del libsbml, scipy, warn
