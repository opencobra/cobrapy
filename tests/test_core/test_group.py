"""Test functions of group via model.py."""

from os.path import join

import pytest

from cobra import Model
from cobra.core import Group
from cobra.io import load_json_model, read_sbml_model, save_json_model


def test_group_add_elements(model: Model) -> None:
    """Test adding elements to a group."""
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

    # the order of members should be the same as the loaded one
    for i in range(num_total_members):
        assert group.members[i] == model.reactions[i]


def test_group_kind() -> None:
    """Test SBML compliance and group kind."""
    group = Group("arbitrary_group1")
    with pytest.raises(ValueError) as excinfo:
        group.kind = "non-SBML compliant group kind"
    assert "Kind can only by one of:" in str(excinfo.value)

    group.kind = "collection"
    assert group.kind == "collection"


def test_read_write_json(data_directory, tmp_path):
    model = read_sbml_model(join(data_directory, "e_coli_core.xml"))
    assert model.groups is not None
    assert len(model.groups) == 10
    assert len(model.groups[0].members) == 6

    path_to_file = join(tmp_path, "group_ecoli.json")
    save_json_model(model, path_to_file)
    model_from_json = load_json_model(path_to_file)
    assert model_from_json.groups is not None
    assert len(model_from_json.groups) == 10
    assert len(model_from_json.groups[0].members) == 6
