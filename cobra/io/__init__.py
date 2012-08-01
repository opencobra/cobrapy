from os import name as __name
try:
    from .sbml import create_cobra_model_from_sbml_file as read_sbml_model
    from .sbml import read_legacy_sbml
    from .sbml import write_cobra_model_to_sbml_file as write_sbml_model
except ImportError, error:
    from warnings import warn
    warn("cobra.io.sbml will not be functional: ImportError %s" % error)


if __name != 'java':
    try:
        from .mat import load_matlab_model
        from .mat import save_matlab_model
    except ImportError, error:
        from warnings import warn
        warn("cobra.io.mat is not be functional: ImportError %s"%error)

del __name



