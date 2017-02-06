try:
    import numpy
except:
    numpy = None

try:
    import scipy
except:
    scipy = None

from .variability import flux_variability_analysis, find_blocked_reactions
from .single_deletion import single_gene_deletion, single_reaction_deletion
from .parsimonious import optimize_minimal_flux
from .loopless import construct_loopless_model
from .gapfilling import growMatch

if numpy:
    from .double_deletion import double_reaction_deletion, double_gene_deletion
    from .phenotype_phase_plane import calculate_phenotype_phase_plane
    from .sampling import sample
else:
    from warnings import warn
    warn("double deletions, phase planes and flux sampling requires numpy")
    del warn
del numpy
