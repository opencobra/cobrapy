import sys
__version__ = '0.1.0'
from os import name as __name
if __name == 'java':
    #raise Exception("Experimental modules of numpy/scipy for java that are" +\
    #"not yet ready for prime time.")
    import oven.danielhyduke.jython.numpy as numpy
    import oven.danielhyduke.jython.scipy as scipy
    from core import Object, Formula, Metabolite, Gene, Reaction, Model
else:
    from core import *
del __name
__doc__ = """
"""
