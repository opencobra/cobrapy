# Release notes for cobrapy x.y.z

## New features

- deletion result DataFrames as returned by functions from `cobra.deletions`
  now have a new `knockout` accessor. See the docs on usage examples.

## Fixes

- remove the frozenset indexing in deletion DataFrames that is now unsupported
  in pandas

## Deprecated features

## Backwards incompatible changes

- deletion result DataFrames have no frozenset index anymore but cow carry the
  deleted elements in the `ids` column.
