# -*- coding: utf-8 -*-

"""Test functions of model.py"""


import pytest

from cobra.core import Group


def test_group_add_elements(model):
    num_members = 5
    reactions_for_group = model.reactions[0:num_members]
    group = Group("arbitrary_group1")
    group.add_members(reactions_for_group)
    group.kind = "collection"
    # number of member sin group should equal the number of reactions
    # assigned to the group
    assert len(group.members) == num_members

    # Choose an overlapping, larger subset of reactions for the group
    num_total_members = 12
    reactions_for_larger_group = model.reactions[0:num_total_members]
    group.add_members(reactions_for_larger_group)
    assert len(group.members) == num_total_members


def test_group_kind():
    group = Group("arbitrary_group1")
    with pytest.raises(ValueError) as excinfo:
        group.kind = "non-SBML compliant group kind"
    assert "Kind can only by one of:" in str(excinfo.value)

    group.kind = "collection"
    assert group.kind == "collection"
