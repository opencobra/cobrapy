# Release notes for cobrapy x.y.z

## New features

- `flux_variability_analysis` gained a `return_fluxes` argument. FVA solves one optimization per reaction per direction and keeps only the objective value; setting this to `True` also returns the flux distribution found at each optimum, as a dict mapping `"minimum"`/`"maximum"` to data frames indexed by the optimized reaction. This avoids re-running an FVA-sized batch of solves when the distributions themselves are wanted, for example as a starting pool for sampling. With `loopless="fastSNP"` the loop constraints are applied to the whole model when this is set, so every returned distribution is loopless. The default return value is unchanged.

## Fixes

## Other

## Deprecated features

## Backwards incompatible changes
