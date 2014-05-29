try:
    import numpy
except:
    numpy = None

from .essentiality import assess_medium_component_essentiality
from .variability import flux_variability_analysis
from .single_deletion import single_deletion

if numpy:
    from .double_deletion import double_deletion
else:
    from warnings import warn
    warn("double_deletion requires numpy")
    del warn
del numpy
