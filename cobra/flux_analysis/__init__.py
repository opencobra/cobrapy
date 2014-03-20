#from .essentiality import assess_medium_component_essentiality
#from .variability import flux_variability_analysis
#from .single_deletion import single_deletion

from os import name as __name
from warnings import warn
if __name == 'java':
    warn('double_deletion functions and moma are not yet supported on %s'%__name)
else:
    try:
        from .double_deletion import double_deletion
    except Exception, e:
        from warnings import warn
        warn("double_deletion is not accessible: %s"%e)