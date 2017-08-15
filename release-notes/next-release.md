# Release notes for cobrapy x.y.z

## Fixes

- `Model.compartment` is now a dynamic property fetching the
  compartments from all metabolites therefore always
  up-to-date. Assigning a dictionary to the same property updates the
  internal dictionary of compartment descriptions. This change removes
  the need for the check for missing compartments from
  `validation.check_metabolite_compartment_formula`.

## New features

## Deprecated features

- `Model.get_metabolite_compartments` is deprecated (use
  `Model.compartments` instead).
- `Reaction.get_compartments` is deprecated (use
  `Reaction.compartments` instead).
