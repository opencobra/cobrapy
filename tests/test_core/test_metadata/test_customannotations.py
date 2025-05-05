"""Test functions of custom.py."""

from pathlib import Path

import pytest

from cobra.core import Model, Species
from cobra.core.metadata.custom import CustomAnnotation, CustomAnnotationStore
from cobra.core.metadata.resource import Qualifier
from cobra.io import (
    load_json_model,
    load_yaml_model,
    read_sbml_model,
    save_json_model,
    save_yaml_model,
    write_sbml_model,
)


def test_customannotation():
    """Test creating and manipulating a single CustomAnnotation."""
    ann = CustomAnnotation.from_data(
        {
            "id": "KV_id",
            "name": "abc_xyz",
            "key": "keyX",
            "value": "45",
            "uri": "https://cobrapy.readthedocs.io/",
        }
    )
    assert isinstance(ann, CustomAnnotation)
    assert ann.id == "KV_id"
    assert ann.name == "abc_xyz"
    assert ann.key == "keyX"
    assert ann.value == "45"
    assert ann.uri == "https://cobrapy.readthedocs.io/"

    with pytest.raises(ValueError):
        ann.remove_from_parent()

    with pytest.raises(AttributeError):
        ann.key = "new_key"
    ann.value = "42"
    assert ann.value == "42"
    ann.uri = "https://opencobra.github.io/"
    assert ann.uri == "https://opencobra.github.io/"

    with pytest.raises(TypeError):
        ann = CustomAnnotation.from_data("key:value")
    ann = CustomAnnotation.from_data(
        {
            "key": "keyX",
            "value": "45",
            "uri": "https://cobrapy.readthedocs.io/",
        }
    )
    assert isinstance(ann.__str__(), str)
    assert isinstance(ann.__repr__(), str)


def test_customannotationstore():
    """Test creating and manipulating a CustomAnnotationStore object."""
    entry1 = {
        "key": "key1",
        "value": "45",
        "uri": "https://cobrapy.readthedocs.io/",
    }
    entry2 = CustomAnnotation.from_data(
        {
            "key": "key2",
            "value": "48",
            "uri": "https://tinyurl2.com/ybyr7b62",
        }
    )
    entry3 = CustomAnnotation(key="key3", value="50")

    kvp = CustomAnnotationStore(entries=[entry1, entry2, entry3])

    assert len(kvp) == 3
    for key in ["key1", "key2", "key3"]:
        assert key in kvp
    assert kvp["key2"] == entry2

    kvp["key1"] = {"key": "key1", "value": "test"}
    assert kvp["key1"].value == "test"
    with pytest.raises(ValueError):
        kvp["key1"] = {"key": "key4", "value": "test"}

    kvp2 = CustomAnnotation(key="key2", value="test")
    kvp["key2"] = kvp2
    assert kvp["key2"].value == "test"
    with pytest.raises(ValueError):
        kvp["key2"] = CustomAnnotation(key="key4", value="test")

    with pytest.raises(ValueError):
        other_store = CustomAnnotationStore(entries=[kvp2])
    other_store = CustomAnnotationStore()
    with pytest.raises(ValueError):
        other_store.add(kvp2)
    with pytest.raises(ValueError):
        other_store[kvp2.key] = kvp2
    with pytest.raises(ValueError):
        kvp[kvp2.key] = kvp2

    kvp["key3"] = "test"
    assert kvp["key3"].value == "test"
    with pytest.raises(TypeError):
        kvp["key1"] = 10

    kvp.add(dict(key="key4", value="value4"))
    assert kvp["key4"].value == "value4"

    with pytest.raises(IndexError):
        kvp.add([dict(key="key4", value="other_value")])

    assert len(kvp) == 4
    kvp.remove(kvp2)
    assert len(kvp) == 3
    with pytest.raises(ValueError):
        kvp.remove("key5")
    with pytest.raises(ValueError):
        kvp.remove(CustomAnnotation(key="key1", value="different"))
    with pytest.raises(ValueError):
        kvp.remove(CustomAnnotation(key="key5", value="different"))

    kvp.remove("key1")
    assert len(kvp) == 2

    kvp["key3"].remove_from_parent()

    assert isinstance(kvp.__repr__(), str)
    assert isinstance(kvp._repr_html_(), str)


def test_customannotation_for_object() -> None:
    """Test accessing custom annotations of a cobrapy object."""
    s = Species()
    s.metadata.custom["key1"] = "value1"
    assert s.metadata.custom["key1"].value == "value1"
    s.add_annotations(
        CustomAnnotation(
            key="key2", value="value2", uri="https://cobrapy.readthedocs.io/"
        )
    )
    assert s.metadata.custom["key2"].value == "value2"
    assert s.metadata.custom["key2"].uri == "https://cobrapy.readthedocs.io/"

    assert s.metadata.to_dict()["custom"]["key2"]["value"] == "value2"


@pytest.mark.parametrize("model_type", ["direct", "sbml", "json", "yaml"])
def test_read_write_model(complex_model: Model, tmp_path: Path, model_type: str):
    """Test annotation consistency when writing and reading an SBML file."""
    model = None
    if model_type == "sbml":
        out_path = tmp_path / "complex_model_writing.sbml"
        assert write_sbml_model(complex_model, str(out_path)) is None

        model = read_sbml_model(str(out_path))
    elif model_type == "json":
        out_path = tmp_path / "complex_model_writing.json"
        assert save_json_model(complex_model, str(out_path)) is None

        model = load_json_model(str(out_path))
    elif model_type == "yaml":
        out_path = tmp_path / "complex_model_writing.yaml"
        assert save_yaml_model(complex_model, str(out_path)) is None

        model = load_yaml_model(str(out_path))
    else:
        model = complex_model
    assert model

    assert (
        model.metadata.standardized[Qualifier.Biological_hasTaxon][0]
        .resources[0]
        .identifier
        == "511145"
    )

    assert model.metadata.custom
    assert model.metadata.custom["trusted_source"].value == "42"

    metabolite = model.metabolites.get_by_id("2pg_c")
    assert metabolite.metadata.custom
    assert metabolite.metadata.custom["toolset"].value == "expert"
    assert (
        metabolite.metadata.custom["toolset"].uri == "urn:awesometool.com:keyvaluepair"
    )

    assert metabolite.metadata.custom["nested"].value == "yes"

    assert (
        metabolite.metadata.custom["nested"]
        .metadata.standardized[Qualifier.Modelling_isDescribedBy][0]
        .resources[0]
        .identifier
        == "e_coli_core"
    )


@pytest.mark.parametrize(
    "model_type",
    [
        "direct",
        pytest.param(
            "sbml", marks=pytest.mark.xfail(reason="not implemented in libsbml")
        ),
        "json",
        "yaml",
    ],
)
def test_read_write_nested_kvp_model(
    complex_model: Model, tmp_path: Path, model_type: str
):
    """Test nested kvp annotation consistency when writing and reading an SBML file."""
    model = None
    if model_type == "sbml":
        out_path = tmp_path / "complex_model_writing.sbml"
        assert write_sbml_model(complex_model, str(out_path)) is None

        model = read_sbml_model(str(out_path))
    elif model_type == "json":
        out_path = tmp_path / "complex_model_writing.json"
        assert save_json_model(complex_model, str(out_path)) is None

        model = load_json_model(str(out_path))
    elif model_type == "yaml":
        out_path = tmp_path / "complex_model_writing.yaml"
        assert save_yaml_model(complex_model, str(out_path)) is None

        model = load_yaml_model(str(out_path))
    else:
        model = complex_model
    assert model

    metabolite = model.metabolites.get_by_id("2pg_c")

    assert (
        metabolite.metadata.custom["nested"].metadata.custom["nested_kvp"].value
        == "not_recommended"
    )
