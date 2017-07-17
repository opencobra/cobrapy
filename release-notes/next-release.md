# Release notes for cobrapy x.y.z

## Fixes

## New features

- `Model.slim_optimize()` can be used perform optimization without
  creating a solution. Can lead to significant speedup compared to
  `Model.optimize` when repeatedly doing optimizations and only making
  use of the objective value as avoiding the need to fetch all values
  from the solver object.
- solution, model, metabolite and reaction now have html
  representation so they give more informative prints in jupyter
  notebooks.
- `Model.medium` now returns a data frame instead of dictionary.

## Deprecated features

