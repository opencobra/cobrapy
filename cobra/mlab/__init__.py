from os import name as __name
from sys import modules as __modules
from warnings import warn
if __name == 'java':
    warn("%s is not yet supported on jython"%__modules[__name__])
else:
    from mlab import *
del __name, __modules

__doc__ = """
Provides python functions which are useful for interacting with MATLAB
by using the mlabwrap library.

To use MATLAB functions directly, see cobra.matlab
"""
