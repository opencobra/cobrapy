# Release notes for cobrapy 0.33.0

## New features

- `spectra_cc` implements the SPECTRA consistency check, a sibling of `fastcc` for removing the blocked reactions of a model. It alternates a forward and a reverse LP, each driving as many of the still-unexplained reactions as possible towards a flux threshold at once, until that set stops shrinking. Unlike FASTCC it needs no reaction flipping: the forward LP's auxiliary variables are unbounded from below, so reactions that can only carry negative flux are simply picked up by the reverse LP instead. The two return the same consistent set on every model tested, and `spectra_cc` is faster than `fastcc` on larger, more reversibility-rich models.

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
