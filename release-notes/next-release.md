# Release notes for cobrapy x.y.z

## New features

## Fixes

Fixed an issue where reaction bounds were being overwritten by global model min/max values
when writing sbml model to file. See [#1300](https://github.com/opencobra/cobrapy/pull/1312).

Fix an issue where [`runfrog`](https://github.com/matthiaskoenig/fbc_curation/issues/98) does
not work via github actions or local installation by removing the use of obsolete numpy
aliases for `float` and `bool`.

## Other

## Deprecated features

## Backwards incompatible changes
