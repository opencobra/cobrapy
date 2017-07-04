# Release notes for cobrapy x.y.z

## Fixes

- `cobra.flux_analysis.reaction.assess`
  [was broken](https://github.com/opencobra/cobrapy/issues/537)
  following the release of 0.6.0 and has now been fixed (and now with
  unit tests).

## New features
- `flux_variability_analysis` now has the `pfba_factor` parameter
  which enables the inclusion of a constraint on the max sum of
  absolute fluxes when doing FVA.

## Deprecated features

- `cobra.flux_analysis.reaction.assess_{precursors,products}` were
  essentially copies of each other and have been merged to
  `cobra.flux_analysis.reaction.assess_component`
