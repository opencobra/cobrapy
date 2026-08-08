# Release notes for cobrapy x.y.z

## New features

- `flux_variability_analysis` gained a `return_fluxes` argument. FVA solves one optimization per reaction per direction and keeps only the objective value; setting this to `True` also returns the flux distribution found at each optimum, as a dict mapping `"minimum"`/`"maximum"` to data frames indexed by the optimized reaction. This avoids re-running an FVA-sized batch of solves when the distributions themselves are wanted, for example as a starting pool for sampling. With `loopless="fastSNP"` the loop constraints are applied to the whole model when this is set, so every returned distribution is loopless. The default return value is unchanged.
- Added `chrr` sampler for flux polytope sampling, which is guaranteed to yield uniformly distributed samples. Uses the optional dependency `hopsy`, installable with `pip install cobra[chrr]`.
- `flux_variability_analysis` gained `time_limit` (a per-solve solver time limit in seconds) and `accept_incumbent`. With `accept_incumbent=True`, a solve interrupted at the limit records the solver's best feasible solution -- flagged per direction in new boolean result columns `minimum_proven`/`maximum_proven` -- instead of NaN; unproven values are one-sided estimates. This makes loopless FVA with the MILP variants usable at genome scale, where a minority of solves cannot prove optimality in any reasonable time: on Recon2 with `loopless="potentials"` and a 4-second limit, roughly half the solves prove optimality and the rest return feasible loopless flux distributions (with `return_fluxes`) suitable for seeding samplers. Non-optimal solves without an incumbent now safely record NaN instead of raising from attribute extraction.

## Fixes
- Rare race condition in cache directory creation from running seperate processes loading cobrapy on clean machine fixed. (https://github.com/opencobra/cobrapy/issues/1476)
- `Formula.__add__` now accepts a plain `str` for the right-hand operand, matching its documented type hint and docstring. Previously `Formula("H2O") + "C"` raised `AttributeError`.


## Other

## Deprecated features


## Backwards incompatible changes
