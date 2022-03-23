# Release notes for cobrapy x.y.z

## New features

Added knock_out_model_genes to simulate knock-outs, setting
reaction bounds of effected reactions to zero and returning a list
of reactions that were knocked-out. Replaces delete_model_genes and undelete_model_genes
since it is context sensitive.
Improve reading of matlab models, which will include metaboilte
and reaction annotations.

## Fixes

`model.copy()` will now correctly copy GPRs.

Fixed an error where matlab models can not be read if their bounds exceed the configuration
default in some cases.

Fixed some bugs in GPR().from_string() where it was using the unmodified string, 
leading to errors with GPRs that should work. Made GPRs that have empty parenthesis 
fail more comprehensibly.

## Other

Moved tests to a different directory, from `src/cobra/test` 
to `src/tests`

Added two tests for GPR fixes
test_gpr_wrong_input()
test_gpr_that_needs_two_replacements()

## Deprecated features

Deprecated delete_model_genes, undelete_model_genes.

## Backwards incompatible changes

removed find_gene_knockout_reactions from delete.py

removed _find_gene_knockout_reactions_fast, 
_gene_knockout_computation, _get_removed 
from test_delete.py
Removed pymatlib direct transfer of models to matlab process. 
Please use save_matlab_model() and then read the model in matlab.

