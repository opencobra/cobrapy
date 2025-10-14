# Release notes for cobrapy x.y.z

## New features

## Fixes

Fixes failures of GPR.copy() in Python 3.13.

Fix compartment not being stored for metabolites created during
reaction.build_reaction_from_string

Fix `reaction.check_mass_balance` giving incorrect results for reactions with floating point coefficients.

Fixes FastCC. This now implements the full algorithm as in the original paper and gives
the same results as the Matlab implementation (within the solver tolerance).

## Other

## Deprecated features

## Backwards incompatible changes

Following libSBML we now also dropped support for Python 3.8. You can still use cobrapy
with Python 3.8 by installing version 0.29.1 or earlier.
