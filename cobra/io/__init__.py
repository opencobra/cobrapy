from warnings import warn

from .sbml3 import read_sbml_model, write_sbml_model
from .json import load_json_model, save_json_model, to_json

# These functions have other dependencies
try:
    import libsbml
except ImportError:
    warn("cobra.io.sbml requires libsbml")
    libsbml = None
else:
    from .sbml import read_legacy_sbml
    from .sbml import write_cobra_model_to_sbml_file as write_legacy_sbml

try:
    import scipy
except ImportError:
    warn("cobra.io.mat requires scipy")
    scipy = None
else:
    from .mat import load_matlab_model, save_matlab_model

del libsbml, scipy, warn
