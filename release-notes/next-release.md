# Release notes for cobrapy x.y.z

## New features

## Fixes

* Reformatted example files (e_coli_core.xml, mini_cobra.xml, mini.json, mini.yml, textbook.xml.gz) to be more compliant with identifiers.org. ncbigi is not a valid identifiers.org, so it was replaced with ncbiprotein.
* make sbml.py subsystem reading add partonomy, which matches the definition
of SBO:0000633 (see https://sourceforge.net/p/sbo/term-request/113/)
* Correct reading and writing of subsystem in mat.
* General cleanup of code in mat.py

## Other

* Resolve `flake8` issues and add missing type annotations and docstrings in `src/cobra/io` and `tests/test_io` (#1212).

## Deprecated features

## Backwards incompatible changes
