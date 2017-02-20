# -*- coding: utf-8 -*-
try:
    import numpy
except:
    numpy = None

try:
    import scipy
except:
    scipy = None

from cobra.flux_analysis.variability import flux_variability_analysis, find_blocked_reactions
from cobra.flux_analysis.single_deletion import single_gene_deletion, single_reaction_deletion
from cobra.flux_analysis.parsimonious import optimize_minimal_flux
from cobra.flux_analysis.loopless import construct_loopless_model
from cobra.flux_analysis.gapfilling import growMatch

if numpy:
    from cobra.flux_analysis.double_deletion import double_reaction_deletion, double_gene_deletion
    from cobra.flux_analysis.phenotype_phase_plane import calculate_phenotype_phase_plane
    from cobra.flux_analysis.sampling import sample
else:
    from warnings import warn
    warn("double deletions, phase planes and flux sampling requires numpy")
    del warn
del numpy
