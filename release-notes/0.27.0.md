# Release notes for cobrapy 0.27.0

## New features

## Fixes

COBRApy is now compatible with pydantic version 2.

COBRApy is now compatible with the latest numpy versions and pin to a lower version
has been removed.

`loopless_solution` now fixes the objective to its optimum as in the
originally published method and returns the objective value in the solution object.

Repair a broken test for `fix_objective_as_constraint`.

## Other

Backwards compatibility for pickled models has been improved.

## Deprecated features

## Backwards incompatible changes
