# Release notes for cobrapy x.y.z

## New features

* Added a new "hybrid" solver that exposes an HIGHS/OSQP combinations for large scale
  LPs, MILPs, and QPs.

## Fixes

Repaired the broken Biomodels Web IO.

## Other

## Deprecated features

* The OSQP solver is deprecated in favor of the hybrid solver which also uses OSQP but
  will not attempt to solve LPs with OSQP anymore. Setting the `model.solver = "osqp"`
  will now use the hybrid interface and will raise an error in a future release.

## Backwards incompatible changes
