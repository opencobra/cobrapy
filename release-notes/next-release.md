# Release notes for cobrapy x.y.z

## New features

Adds a `fast_gapfill` parameter to `cobra.flux_analysis.gapfilling.gapfill()` to enable Fast Gap Filling, which is significantly faster than the default method, at the cost of some accuracy. This is based on the method described in the supplementary information of the following paper:

> Dreyfuss, J. M., Zucker, J. D., Hood, H. M., Ocasio, L. R., Sachs, M. S., & Galagan, J. E. (2013). [Reconstruction and validation of a genome-scale metabolic model for the filamentous fungus Neurospora crassa using FARM](https://doi.org/10.1371/journal.pcbi.1003126). PLoS Computational Biology, 9(7), e1003126. https://doi.org/10.1371/journal.pcbi.1003126


## Fixes

Fixes failures of GPR.copy() in Python 3.13.


Fix compartment not being stored for metabolites created during
reaction.build_reaction_from_string

Fix `reaction.check_mass_balance` giving incorrect results for reactions with floating point coefficients.

## Other

## Deprecated features

## Backwards incompatible changes

Following libSBML we now also dropped support for Python 3.8. You can still use cobrapy
with Python 3.8 by installing version 0.29.1 or earlier.
