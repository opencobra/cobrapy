try:
    import numpy
except:
    numpy = None

from .essentiality import assess_medium_component_essentiality
from .variability import flux_variability_analysis, find_blocked_reactions
from .single_deletion import single_gene_deletion, single_reaction_deletion
from .parsimonious import optimize_minimal_flux
from .loopless import construct_loopless_model

if numpy:
    from .double_deletion import double_reaction_deletion, double_gene_deletion
else:
    from warnings import warn
    warn("double_deletion requires numpy")
    del warn
del numpy
