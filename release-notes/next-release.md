# Release notes for cobrapy x.y.z

## New features

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
