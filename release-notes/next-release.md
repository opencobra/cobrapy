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
- New convenience functions `cobra.flux_analysis.find_essential_genes` and
  `cobra.flux_analysis.find_essential_reactions`.
- `Model.optimize` has new parameter `raise_error` to enable option to
  get trigger exception if no feasible solution could be found.
- `str(reaction)` now gives the more useful reaction id and the
  reaction string.

## Deprecated features

- `str(reaction)` no longer gives `reaction.id`.
