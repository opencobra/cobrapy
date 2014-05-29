from warnings import warn

try:
    import libsbml
except ImportError:
    libsbml = None
try:
    import scipy
except ImportError:
    scipy = None

if libsbml:
    from .sbml import create_cobra_model_from_sbml_file as read_sbml_model
    from .sbml import read_legacy_sbml
    from .sbml import write_cobra_model_to_sbml_file as write_sbml_model
else:
    warn("cobra.io.sbml requires libsbml")

if scipy:
    from .mat import load_matlab_model
    from .mat import save_matlab_model
else:
    warn("cobra.io.mat requires scipy")


from .json import load_json_model
from .json import save_json_model, to_json

del libsbml, scipy, warn
