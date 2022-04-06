# Release notes for cobrapy x.y.z

## New features

Improve reading of matlab models, which will include metaboilte
and reaction annotations.

## Fixes

`model.copy()` will now correctly copy GPRs.

Fixed an error where matlab models can not be read if their bounds exceed the configuration
default in some cases.

## Other

Removed pymatlib direct transfer of models to matlab process. 
Please use save_matlab_model() and then read the model in matlab.

## Deprecated features

## Backwards incompatible changes
