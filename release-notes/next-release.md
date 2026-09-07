# Release notes for cobrapy 0.32.1

## New features

- `flux_variability_analysis` gained an `all_fluxes` argument. FVA solves one optimization per reaction per direction and keeps only the objective value; setting this to `True` also returns the flux distribution found at each optimum, as a dict mapping `"minimum"`/`"maximum"` to data frames indexed by the optimized reaction. This avoids re-running an FVA-sized batch of solves when the distributions themselves are wanted, for example as a starting pool for sampling. With `loopless="fastSNP"` the loop constraints are applied to the whole model when this is set, so every returned distribution is loopless. The default return value is unchanged.

- `Diitlists` (`model.reactions`, `model.genes`, `model.metabolites`) are now generic containers and the types of their contents are now inferable by python type checkers.

## Fixes

- fixed import failure when hopsy was not installed

## Other

## Deprecated features

## Backwards incompatible changes
