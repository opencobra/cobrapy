# Release notes for cobrapy x.y.z

## New features

## Fixes

* serializes GPRs to strings to avoid massive storage usage
* Reformatted example files (e_coli_core.xml, mini_cobra.xml, mini.json, mini.yml, textbook.xml.gz) to be more compliant with identifiers.org. ncbigi is not a valid identifiers.org, so it was replaced with ncbiprotein.
* make sbml.py subsystem reading add partonomy, which matches the definition
of SBO:0000633 (see https://sourceforge.net/p/sbo/term-request/113/)
* Correct reading and writing of subsystem in mat.
* General cleanup of code in mat.py
* Fixed issue 673 (problem when copying reaction from a different model)

## Other

* Resolve `flake8` issues and add missing type annotations and docstrings in `src/cobra/io` and `tests/test_io` (#1212).
* Updated model.py and test_model.py to Python 3.6+, including type annotations and docstrings.

## Deprecated features

## Backwards incompatible changes
* Removed `model.add_reaction()` and replaced remaining usages of it with `model.add_reactions()`
* Removed the following tests: test_add_remove_reaction_benchmark, test_add_reaction, test_add_reaction_context, test_add_reaction_from_other_model, test_add_cobra_reaction
* Removed `model.__add__` and `model.__iadd__` - use `model.merge` to replace them.
* Remove `Model().description()`.
* Remove `Model().get_metabolite_compartments()`.
