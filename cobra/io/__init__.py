from warnings import warn
try:
    from sbml import *
except:
    warn ('Unable to import cobra.io.sbml.  Perhaps your external sbml libraries are not installed?')
