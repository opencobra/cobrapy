# Release notes for cobrapy x.y.z

## New features

- Added `nullspace_fast_snp`, a new function that constructs sparse basis for the nullspace of a matrix using Fast-SNP algorithm.
- Added Fast-SNP based method of adding loopless constraints to `add_loopless` and `flux_variability_analysis`.
- Added `find_cyclic_reactions`, a new function that finds all reactions in a model that can be part of a cycle.
- Optimized `flux_variability_analysis` by running MILP optimization only for reactions that can be a part of cycle (as identified by `find_cyclic_reactions`).

## Fixes

- Updated BioModel URL from `https://www.ebi.ac.ak/biomodels/model/` to `https://biomodels.org/model/`.

## Other

- Updated type hinting for the `DictList` class so that the type of `Object` contained by a `DictList` can be specified. For example, the hinted return type of `model.reactions.get_by_id` is now `Reaction` instead of `Object`.

## Deprecated features

- Changed the type of the `loopless` parameter in `flux_variability_analysis` from `bool` to `Optional[str]`. Using `loopless=False` or `loopless=True` (boolean) is now deprecated.

## Backwards incompatible changes
