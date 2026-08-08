# Release notes for cobrapy x.y.z

## New features

- `add_loopless` gained `method="potentials"`, which encodes the loop law through free metabolite potentials (`G = S_intᵀμ`) instead of a computed null-space basis. The feasible flux space is identical to the existing methods; only construction cost differs. On Recon2 the null-space step dominates construction (98% of a 17-minute build), which "potentials" skips entirely, building the same constraints in seconds. `flux_variability_analysis` accepts `loopless="potentials"` accordingly.

## Fixes

## Other

## Deprecated features

## Backwards incompatible changes
