import sys
__version__ = '0.1.0.b1'
if hasattr(sys, 'JYTHON_JAR'):
    raise Exception("Experimental modules of numpy/scipy for java that are" +\
    "not yet ready for prime time.")
    import cobra.oven.hyduke.numpy as numpy
    import cobra.oven.hyduke.scipy as scipy
from core import *
__doc__ = """
"""
