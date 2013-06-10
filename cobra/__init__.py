import sys
__version__ = '0.3.0-dev'
from os import name as __name
from core import Object, Formula, Metabolite, Gene, Reaction, Model, DictList, Species
if __name == 'java':
    #raise Exception("Experimental modules of numpy/scipy for java that are" +\
    #"not yet ready for prime time.")
    #import oven.danielhyduke.jython.numpy as numpy
    #import oven.danielhyduke.jython.scipy as scipy
    from warnings import warn
    warn("COBRA for Python is not optimized for JAVA.  If it's slow or crashes consider increasing JVM memory")
else:
    try:
        from core import ArrayBasedModel
    except Exception, e:
        from warnings import warn
        warn("cobra.ArrayBasedModel class is unavailable: %s"%repr(e))

del __name
__doc__ = """
"""
