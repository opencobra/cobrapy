try:
    from essentiality import assess_medium_component_essentiality
except Exception, e:
    from warnings import warn
    warn("essentiality is not accessible: %s"%e)
try:
    from variability import flux_variability_analysis
except Exception, e:
    from warnings import warn
    warn("variability is not accessible: %s"%e)
try:
    from single_deletion import single_deletion
except Exception, e:
    from warnings import warn
    warn("single_deletion is not accessible: %s"%e)


from os import name as __name
from warnings import warn
if __name == 'java':
    warn('double_deletion functions and moma are not yet supported on %s'%__name)
else:
    try:
        from double_deletion import double_deletion
    except Exception, e:
        from warnings import warn
        warn("double_deletion is not accessible: %s"%e)

    







