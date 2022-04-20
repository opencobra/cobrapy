# Release notes for cobrapy x.y.z

## New features

## Fixes

`model.copy()` will now correctly copy GPRs.

Fixed an error where matlab models can not be read if their bounds exceed the configuration
default in some cases.

Fixed some bugs in GPR().from_string() where it was using the unmodified string, 
leading to errors with GPRs that should work. Made GPRs that have empty parenthesis 
fail more comprehensibly.

## Other

Added two tests for GPR fixes
test_gpr_wrong_input()
test_gpr_that_needs_two_replacements()

## Deprecated features

## Backwards incompatible changes
