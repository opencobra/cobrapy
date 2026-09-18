# Release notes for cobrapy 0.33.0

## New features

- `pfba` and `geometric_fba` gained a `raise_error` argument that is forwarded to the underlying solver-status check, so all FBA flavours can be invoked like `Model.optimize(raise_error=...)`. As with `Model.optimize`, a non-optimal status now emits a warning by default and only raises when `raise_error=True` (previously `pfba` always raised).

- `flux_variability_analysis` gained an `all_fluxes` argument. FVA solves one optimization per reaction per direction and keeps only the objective value; setting this to `True` also returns the flux distribution found at each optimum, as a dict mapping `"minimum"`/`"maximum"` to data frames indexed by the optimized reaction. This avoids re-running an FVA-sized batch of solves when the distributions themselves are wanted, for example as a starting pool for sampling. With `loopless="fastSNP"` the loop constraints are applied to the whole model when this is set, so every returned distribution is loopless. The default return value is unchanged.

- `Dictlists` (`model.reactions`, `model.genes`, `model.metabolites`) are now generic containers and the types of their contents are now inferable by python type checkers.

## Fixes

- type-annotations on resetable properties fixed.
- fixed the URL for the BIGG repository
- fixed the URL for the BioModels repository

## Other

- Downloads from BioModels and BIGG now get xfail markers as they are often blocked

## Deprecated features

## Backwards incompatible changes
