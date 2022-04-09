# Release notes for cobrapy x.y.z

## New features

Added knock_out_model_genes to simulate knock-outs, setting
reaction bounds of effected reactions to zero and returning a list
of reactions that were knocked-out. Replaces delete_model_genes and undelete_model_genes
since it is context sensitive.

## Fixes

`model.copy()` will now correctly copy GPRs.

Fixed an error where matlab models can not be read if their bounds exceed the configuration
default in some cases.

## Other

## Deprecated features

Deprecated delete_model_genes, undelete_model_genes.

## Backwards incompatible changes

removed find_gene_knockout_reactions from delete.py

removed _find_gene_knockout_reactions_fast, 
_gene_knockout_computation, _get_removed 
from test_delete.py
