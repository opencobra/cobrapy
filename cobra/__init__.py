import sys
from .version import get_version
__version__ = get_version()
from os import name as __name
from .core import Object, Formula, Metabolite, Gene, Reaction, Model, DictList, Species
if __name == 'java':
    from warnings import warn
    warn("COBRA for Python is not optimized for JAVA. If it's slow or crashes consider increasing JVM memory")
else:
    try:
        from .core import ArrayBasedModel
    except Exception as e:
        None

del __name, get_version
__doc__ = """
"""
