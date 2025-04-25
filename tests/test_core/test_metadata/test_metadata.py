"""Tests for the metadata structures."""

import json
from pathlib import Path

import pytest

from cobra import Model
from cobra.core.metadata import Metadata, Qualifier, StandardizedAnnotation
from cobra.core.metadata.resource import QualifiersAlias, Resource
from cobra.core.metadata.standardized import (
    SimplifiedAnnotationInterface,
    StandardizedAnnotationList,
    StandardizedAnnotationStore,
)
from cobra.core.species import Species
from cobra.io import load_json_model, read_sbml_model, save_json_model, write_sbml_model


PUBMED_EXAMPLE = "https://identifiers.org/pubmed/1111111"
ECO_EXAMPLE = "https://identifiers.org/eco/ECO:0000004"
RESOURCE_LIST = [
    "http://identifiers.org/bigg.metabolite/13dpg",
    "http://identifiers.org/biocyc/DPG",
    "http://identifiers.org/chebi/CHEBI:11881",
    "http://identifiers.org/chebi/CHEBI:16001",
    "http://identifiers.org/chebi/CHEBI:1658",
    "http://identifiers.org/chebi/CHEBI:20189",
    "http://identifiers.org/chebi/CHEBI:57604",
    "http://identifiers.org/hmdb/HMDB01270",
    "http://identifiers.org/kegg.compound/C00236",
    "http://identifiers.org/pubchem.substance/3535",
    "http://identifiers.org/reactome/REACT_29800",
    "http://identifiers.org/seed.compound/cpd00203",
    "http://identifiers.org/unipathway.compound/UPC00236",
]
ECOLI_MODEL_ANNOTATIONS = [
    {
        "qualifier": "bqb_hasTaxon",
        "resources": ["http://identifiers.org/taxonomy/511145"],
    },
    {
        "qualifier": "bqm_is",
        "resources": ["http://identifiers.org/bigg.model/e_coli_core"],
        "annotations": [
            {
                "qualifier": "bqb_isDescribedBy",
                "resources": [PUBMED_EXAMPLE],
            },
            {
                "qualifier": "bqb_isDescribedBy",
                "resources": [ECO_EXAMPLE],
            },
        ],
    },
    {
        "qualifier": "bqm_isDescribedBy",
        "resources": ["http://identifiers.org/doi/10.1128/ecosalplus.10.2.1"],
    },
    {
        "qualifier": "bqm_isDescribedBy",
        "resources": ["http://identifiers.org/ncbiprotein/16128336"],
    },
]
COBRA_URL = "https://cobrapy.readthedocs.io/"
CHEBI_SET = {"CHEBI:43215", "CHEBI:11881"}


def test_annotation() -> None:
    """Test creating an annotation manually."""
    # a cobra component
    s = Species()
    # assert s.metadata == {}  # nothing set for annotation, so empty dict
    assert s.metadata.standardized == StandardizedAnnotationStore()
    # assert not s.annotation.keys()
    # assert s.metadata.custom == {}
    assert s.metadata.history.creators == []
    assert s.metadata.history.modified_dates == []

    # setting annotation via old annotation format
    s.annotation["chebi"] = list(CHEBI_SET)

    assert set(s.annotation["chebi"]) == CHEBI_SET
    assert set(s.annotation["chebi"]) == CHEBI_SET
    # assert set(s.metadata.standardized[Qualifier.Biological_is]["chebi"]) == CHEBI_SET

    s.metadata.standardized = StandardizedAnnotationStore()

    s.add_annotations(
        [
            "https://identifiers.org/chebi/CHEBI:43215",
            "https://identifiers.org/chebi/CHEBI:11881",
        ]
    )
    assert s.metadata.standardized == {
        "https://identifiers.org/chebi/CHEBI:43215",
        "https://identifiers.org/chebi/CHEBI:11881",
    }

    # checking old (fixed) annotation format
    assert s.annotation == {"chebi": list(sorted(["CHEBI:43215", "CHEBI:11881"]))}

    # checking new standardized
    cvt = StandardizedAnnotationStore(
        [
            StandardizedAnnotation(
                qualifier=Qualifier.Biological_is,
                resources=[
                    "https://identifiers.org/chebi/CHEBI:43215",
                    "https://identifiers.org/chebi/CHEBI:11881",
                ],
            ),
        ]
    )

    # The next assertion should probably not hold, because of how we add annotations
    # with the same qualifier. Or we should change the behaviour of add_annotations
    # etc. when only strings are supplied.
    # assert s.metadata.standardized == cvt
    s.metadata.standardized = []
    assert s.metadata.standardized == StandardizedAnnotationStore()
    assert s.metadata.standardized.resources == frozenset()
    assert s.metadata.standardized == {}

    s.metadata.standardized = cvt

    assert s.metadata.standardized.resources == {
        "https://identifiers.org/chebi/CHEBI:43215",
        "https://identifiers.org/chebi/CHEBI:11881",
    }

    # checking old (fixed) annotation format
    assert s.annotation == {"chebi": list(sorted(["CHEBI:43215", "CHEBI:11881"]))}

    cvt[0].remove_from_parent()
    assert s.metadata.standardized == []

    # adding an SBO term
    s.metadata.sbo = ["SBO:0000123"]
    assert "sbo" in s.annotation
    # assert s.annotation == {
    #     "chebi": ["CHEBI:43215", "CHEBI:11881"],
    #     "sbo": ["SBO:0000123"],
    # }

    cvt2 = StandardizedAnnotationStore(
        [
            StandardizedAnnotation(
                qualifier="bqb_is",
                resources=["https://identifiers.org/chebi/CHEBI:11881"],
            ),
            StandardizedAnnotation("https://identifiers.org/chebi/CHEBI:43215"),
        ]
    )

    s.metadata.standardized = cvt2
    assert s.metadata.standardized.uris == {
        "https://identifiers.org/chebi/CHEBI:43215",
        "https://identifiers.org/chebi/CHEBI:11881",
    }

    # print(s.annotation[Qualifier.Biological_is])
    # print(s.annotation[0])
    print(s.annotation["chebi"])

    assert (
        s.metadata.standardized[0].resources[0].uri
        == "https://identifiers.org/chebi/CHEBI:11881"
    )

    assert set(s.metadata.standardized[Qualifier.Biological_is].uris) == {
        "https://identifiers.org/chebi/CHEBI:43215",
        "https://identifiers.org/chebi/CHEBI:11881",
    }

    # assert 0 == 1
    # s.annotation.__delitem__("sbo")

    # checking old (fixed) annotation format
    # assert s.annotation == {"chebi": sorted(["CHEBI:43215", "CHEBI:11881"])}


def test_standardized_annotation() -> None:
    """Test creating and manipulating standardized annotations."""

    s = Species()
    s.metadata.add_standardized(ECOLI_MODEL_ANNOTATIONS)
    assert s.metadata.standardized == ECOLI_MODEL_ANNOTATIONS

    s_other = Species()
    with pytest.raises(ValueError):
        s_other.add_annotations(s.metadata.standardized[0])

    ann = StandardizedAnnotation("https://identifiers.org/chebi/CHEBI:43215")
    s_other.add_annotations([ann])

    assert (
        s_other.metadata.standardized.resources_for("chebi")[0].identifier
        == "CHEBI:43215"
    )

    ann.remove_from_parent()
    assert len(s_other.metadata.standardized.resources_for("chebi")) == 0
    with pytest.raises(ValueError):
        ann.remove_from_parent()

    ann.qualifier = Qualifier.Biological_hasPart
    assert ann.qualifier == Qualifier.Biological_hasPart
    assert ann.qualifier in QualifiersAlias.Biological_any

    ann.qualifier = "bqm_is"
    assert ann.qualifier == Qualifier.Modelling_is
    assert ann.qualifier not in QualifiersAlias.Biological_known

    ann.resources = [
        "https://identifiers.org/chebi/CHEBI:43215",
        "https://identifiers.org/chebi/CHEBI:11881",
    ]
    assert next(iter(ann.resources)).identifier == "CHEBI:43215"

    ann.annotations = [
        StandardizedAnnotation(["http://identifiers.org/taxonomy/511145"])
    ]
    assert len(ann.uris) == 2
    assert len(ann.all_uris) == 3

    assert all(resource.namespace == "chebi" for resource in ann.resources)
    assert all(namespace == "chebi" for (namespace, _) in ann.to_tuples())

    records = ann.to_records()
    assert len(records) == 3
    main_group = next(x["annotation_group"] for x in records if x["parent_group"] == 0)
    assert (
        next(x["namespace"] for x in records if x["parent_group"] == main_group)
        == "taxonomy"
    )

    ann_dict = ann.to_dict()
    assert ann_dict == ann
    del ann_dict["annotations"]
    assert not ann_dict == ann
    ann_dict = ann.to_dict()
    resource = ann.resources[0]
    with pytest.raises(ValueError):
        ann.add_resources([resource])
    resource.remove_from_parent()
    assert not ann_dict == ann
    with pytest.raises(ValueError):
        resource.remove_from_parent()
    ann.add_resources([resource])
    assert ann_dict == ann
    with pytest.raises(ValueError):
        new_ann = StandardizedAnnotation()
        new_ann.add_resources([resource])

    ann.resources = None
    assert len(ann.resources) == 0

    ann.add_resources(
        [Resource.from_data({"uri": "http://identifiers.org/taxonomy/51114"})]
    )
    assert len(ann.resources) == 1
    ann.add_resources(
        [Resource.from_data({"namespace": "chebi", "identifier": "CHEBI:43215"})]
    )
    assert len(ann.resources) == 2
    with pytest.raises(TypeError):
        ann.add_resources([Resource.from_data(None)])

    resource = Resource(uri="https://identifiers.org/uniprot/A0PK11")

    assert resource.namespace == "uniprot"
    assert (
        Resource(uri="https://www.uniprot.org/uniprotkb/A0PK11", strict=False).namespace
        is None
    )
    with pytest.raises(ValueError):
        Resource(uri="https://www.uniprot.org/uniprotkb/A0PK11")

    assert resource.to_dict()["namespace"] == "uniprot"
    assert resource == Resource(resource.uri)
    assert resource == {"uri": resource.uri}
    assert resource == ("uniprot", "A0PK11")

    assert resource != ("uniprot",)
    assert resource != 1

    with pytest.raises(TypeError):
        ann.resources = [1, 2]

    with pytest.raises(TypeError):
        ann.qualifier = "unknown"

    with pytest.raises(TypeError):
        ann.qualifier = 1

    assert isinstance(ann._repr_html_(), str)
    assert isinstance(resource._repr_html_(), str)


def test_standardized_annotation_store() -> None:
    """Test creating and manipulating standardized stores."""

    s = Species()
    s.metadata.add_standardized(ECOLI_MODEL_ANNOTATIONS)
    assert s.metadata.standardized == ECOLI_MODEL_ANNOTATIONS
    assert isinstance(s.metadata.standardized, StandardizedAnnotationStore)
    s_other = Species()
    s_other.metadata.add_standardized(ECOLI_MODEL_ANNOTATIONS)
    assert s.metadata.standardized == s_other.metadata.standardized
    assert s.metadata == s_other.metadata
    with pytest.raises(TypeError):
        _ = s.metadata == 1

    with pytest.raises(ValueError):
        s.metadata.standardized.add(s.metadata.standardized[0])
    new_store = StandardizedAnnotationStore()
    assert not new_store
    with pytest.raises(ValueError):
        new_store.add(s.metadata.standardized[0])
    new_list = StandardizedAnnotationList()
    assert not new_list
    new_list.add(s.metadata.standardized[0])
    assert len(new_list) == 1

    new_store = StandardizedAnnotationStore(
        [
            StandardizedAnnotation("https://identifiers.org/chebi/CHEBI:43215"),
            "https://identifiers.org/chebi/CHEBI:11881",
            Resource("https://identifiers.org/chebi/CHEBI:16001"),
            {
                "qualifier": Qualifier.Biological_hasTaxon,
                "resources": ["http://identifiers.org/taxonomy/511145"],
            },
        ]
    )
    assert len(new_store) == 3
    assert len(new_store.resources) == 4

    with pytest.raises(TypeError):
        StandardizedAnnotationStore(1)
    with pytest.raises(TypeError):
        StandardizedAnnotationStore([1])

    assert len(s.metadata.standardized.resources) == 4
    assert len(s.metadata.standardized.all_resources) == 6
    assert len(s.metadata.standardized.uris) == 4
    assert len(s.metadata.standardized.all_uris) == 6
    assert len(s.metadata.standardized.qualifiers) == 3
    assert len(s.metadata.standardized.all_qualifiers) == 4

    assert len(s.metadata.standardized[Qualifier.Modelling_isDescribedBy]) == 2
    assert len(s.metadata.standardized[QualifiersAlias.Modelling_any]) == 3
    assert (
        len(
            s.metadata.standardized[
                [Qualifier.Biological_hasTaxon, Qualifier.Modelling_is]
            ]
        )
        == 2
    )
    with pytest.raises(TypeError):
        s.metadata.standardized["a"]

    assert (
        s.metadata.standardized.resources_for(
            qualifier=Qualifier.Modelling_isDescribedBy, namespace="doi"
        )[0].identifier
        == "10.1128/ecosalplus.10.2.1"
    )
    assert (
        s.metadata.standardized.resources_for(
            qualifier=QualifiersAlias.Modelling_known,
            namespace=["eco", "taxonomy"],
            nested=True,
        )[0].uri
        == ECO_EXAMPLE
    )
    with pytest.raises(TypeError):
        assert (
            s.metadata.standardized.resources_for(
                qualifier=[QualifiersAlias.Modelling_known, 1],
                namespace=["eco", "taxonomy"],
                nested=True,
            )[0].uri
            == ECO_EXAMPLE
        )
    ann = new_store[0]
    ann.remove_from_parent()

    s.metadata.standardized[0] = ann
    assert s.metadata.standardized[0] == ann

    assert not (s.metadata.standardized == "string")
    assert s.metadata.standardized != s_other.metadata.standardized
    assert s.metadata != s_other.metadata

    assert isinstance(s.metadata.standardized._repr_html_(), str)


def test_old_style_annotation() -> None:
    """Test creating old style annotations using add_simple_annotations."""
    s = Species()
    s.annotation.add({"chebi": "CHEBI:17234"})
    s.annotation.add({"chebi": ["CHEBI:1723456", "CHEBI:172345"]})
    with pytest.raises(TypeError):
        s.annotation.add({"chebi": [["CHEBI:123", "CHEBI:1234"]]})
    assert len(s.annotation) == 1
    assert s.annotation.number_of_resources == 3
    s.annotation["eco"] = "123"
    assert len(s.annotation) == 2
    assert s.annotation.number_of_resources == 4
    assert set(s.annotation["eco"]) == {"123"}
    assert set(s.annotation.get("eco")) == {"123"}
    assert set(s.annotation.get("eco", ["456"])) == {"123"}
    with pytest.raises(IndexError):
        _ = s.annotation["invalid"]
    assert s.annotation.get("invalid") is None
    assert set(s.annotation.get("invalid", ["456"])) == {"456"}

    s.annotation.update({"bigg.metabolite": "glc__D"})
    assert set(s.annotation["bigg.metabolite"]) == {"glc__D"}
    del s.annotation["bigg.metabolite"]
    assert s.annotation.get("bigg.metabolite") is None
    other_interface = SimplifiedAnnotationInterface(Metadata())
    other_interface["bigg.metabolite"] = "glc__D"
    s.annotation.update(other_interface)
    assert set(s.annotation["bigg.metabolite"]) == {"glc__D"}
    s.annotation.update({"bigg.metabolite": ["glc__L"]})
    assert set(s.annotation["bigg.metabolite"]) == {"glc__L"}
    del s.annotation["bigg.metabolite"]

    ref_ann = Metadata()
    simpl_ann = SimplifiedAnnotationInterface(ref_ann)
    simpl_ann.add(
        {
            "chebi": ["CHEBI:17234", "CHEBI:1723456", "CHEBI:172345"],
            "eco": ["123"],
        }
    )
    assert s.metadata == ref_ann

    s.annotation.delete_annotation("CHEBI:172345")
    assert len(s.annotation) == 2
    assert s.annotation.number_of_resources == 3
    ref_ann = Metadata()
    simpl_ann = SimplifiedAnnotationInterface(ref_ann)
    simpl_ann.add(
        {
            "chebi": ["CHEBI:17234", "CHEBI:1723456"],
            "eco": ["123"],
        }
    )
    assert s.metadata == ref_ann

    print(s.annotation.to_dict())
    del s.annotation["chebi"]
    print(s.annotation.to_dict())
    assert len(s.annotation) == 1
    assert s.annotation.number_of_resources == 1

    s.annotation.add({"chebi": "CHEBI:17234"})
    s.annotation.add({"chebi": ["CHBEI:1723456", "CHEBI:172345"]})
    assert len(s.annotation) == 2
    assert s.annotation.number_of_resources == 4
    s.annotation["chebi"] = ["CHEBI:123", "CHEBI:1234"]
    assert len(s.annotation) == 2
    assert s.annotation.number_of_resources == 3
    ref_ann = Metadata()
    simpl_ann = SimplifiedAnnotationInterface(ref_ann)
    simpl_ann.add({"chebi": ["CHEBI:123", "CHEBI:1234"], "eco": ["123"]})
    assert s.metadata == ref_ann

    assert len(s.annotation.keys()) == 2

    s.annotation["chebi"] = []
    assert s.annotation.number_of_resources == 1

    s.annotation.clear()

    assert len(s.annotation.keys()) == 0

    s.annotation = {"chebi": ["CHEBI:123", "CHEBI:1234"], "eco": ["123"]}
    assert s.annotation.number_of_resources == 3
    s.annotation = {
        "chebi": ["CHEBI:123", "CHEBI:1234"],
        "eco": ["123"],
        "sbo": ["SBO:0000123"],
    }
    assert s.annotation.number_of_resources == 4


def test_nested_annotation(data_directory: Path) -> None:
    """Test reading annotation from JSON, including nested data.

    Parameters
    ----------
    data_directory: Path
    """
    # testing via standardized
    with data_directory.joinpath("cvterms_nested.json").open("r") as f_cvterms:
        cvterms_data = json.load(f_cvterms)

    s = Species()
    s.metadata.add_standardized(cvterms_data)
    # assert s.annotation == {
    #     "chebi": ["CHEBI:17627"],
    #     "eco": ["000000"],
    #     "kegg.compound": ["C00032"],
    #     "pubmed": ["1111111"],
    #     "uniprot": ["P68871", "P69905"],
    # }
    # check standardized
    main_cvt = [
        {
            "resources": [
                "https://identifiers.org/uniprot/P69905",
                "https://identifiers.org/uniprot/P68871",
                "https://identifiers.org/kegg.compound/C00032",
            ],
            "qualifier": "bqb_hasPart",
        },
        {
            "qualifier": "bqb_hasPart",
            "resources": [
                "https://identifiers.org/uniprot/P69905",
                "https://www.uniprot.org/uniprot/P68871",
                "https://identifiers.org/chebi/CHEBI:17627",
            ],
            "annotations": [
                {
                    "qualifier": "bqb_isDescribedBy",
                    "resources": [
                        PUBMED_EXAMPLE,
                        "https://identifiers.org/eco/000000",
                    ],
                }
            ],
        },
    ]
    nested_cvt = [
        {
            "qualifier": "bqb_isDescribedBy",
            "resources": [PUBMED_EXAMPLE, "https://identifiers.org/eco/000000"],
        }
    ]
    assert s.metadata.standardized == main_cvt
    nested_data = s.metadata.standardized[1].annotations
    assert nested_data == nested_cvt

    additional_cvt = {
        "qualifier": "bqm_is",
        "resources": ["https://identifiers.org/bigg.metabolite/hemoglobin"],
    }
    s.metadata.standardized.add(additional_cvt)
    assert (
        len(
            s.metadata.standardized.resources_for(
                qualifier=Qualifier.Biological_hasPart,
            )
        )
        == 6
    )
    assert (
        len(
            s.metadata.standardized.resources_for(
                qualifier=Qualifier.Biological_hasPart, namespace=["uniprot"]
            )
        )
        == 3
    )
    assert (
        len(
            s.metadata.standardized.resources_for(
                qualifier=[Qualifier.Biological_hasPart],
                nested=True,
            )
        )
        == 8
    )
    assert (
        len(
            s.metadata.standardized.resources_for(
                qualifier=Qualifier.Biological_hasPart,
                namespace="pubmed",
                nested=True,
            )
        )
        == 1
    )
    assert (
        len(
            s.metadata.standardized.resources_for(
                qualifier=QualifiersAlias.Any_is,
                nested=True,
            )
        )
        == 1
    )


def test_cvterms_from_ecoli_xml(annotation_model: Model) -> None:
    """Test the new and old style annotations of an ecoli model."""
    qualifier_set = {
        Qualifier(qual) for qual in ["bqb_hasTaxon", "bqm_is", "bqm_isDescribedBy"]
    }
    nested_cvt = [
        {
            "qualifier": "bqb_isDescribedBy",
            "resources": [PUBMED_EXAMPLE],
        },
        {
            "qualifier": "bqb_isDescribedBy",
            "resources": [ECO_EXAMPLE],
        },
    ]
    ecoli_model_cvterm = StandardizedAnnotationStore.from_data(ECOLI_MODEL_ANNOTATIONS)
    print(ecoli_model_cvterm.to_list_of_dicts())
    print(annotation_model.metadata.standardized.to_list_of_dicts())
    xml_model_cvterms = annotation_model.metadata.standardized
    model_cvterms_qualifier_set = xml_model_cvterms.qualifiers
    assert qualifier_set == model_cvterms_qualifier_set
    assert xml_model_cvterms == ecoli_model_cvterm
    assert (
        len(
            annotation_model.metadata.standardized.query(
                "bqm_isDescribedBy", "qualifier"
            )
        )
        == 2
    )
    nested_data = annotation_model.metadata.standardized.query("bqm_is", "qualifier")[
        0
    ].annotations
    assert nested_data == nested_cvt

    # check backwards compatibility
    # assert annotation_model.annotation.metadata == {
    #     "bigg.model": ["e_coli_core"],
    #     "doi": ["10.1128/ecosalplus.10.2.1"],
    #     "eco": ["ECO:0000004"],
    #     "ncbiprotein": ["16128336"],
    #     "pubmed": ["1111111"],
    #     "taxonomy": ["511145"],
    # }
    # annotation_model.annotation.standardized.delete_annotation("coli")
    # assert annotation_model.annotation.metadata == {
    #     "doi": ["10.1128/ecosalplus.10.2.1"],
    #     "eco": ["ECO:0000004"],
    #     "ncbiprotein": ["16128336"],
    #     "pubmed": ["1111111"],
    #     "taxonomy": ["511145"],
    # }


def test_writing_xml(annotation_model: Model, tmp_path):
    """Test writing a model with annotations to xml (SBML)."""
    assert (
        write_sbml_model(
            annotation_model, str(tmp_path.joinpath("e_coli_core_writing.xml"))
        )
        is None
    )
    # TODO: Add more tests here.


def test_read_write_json(annotation_model: Model, tmp_path: Path):
    """Test writing a model with annotations to JSON."""
    json_path = tmp_path / "e_coli_core_json_writing.json"
    print(json_path)
    assert save_json_model(annotation_model, json_path, sort=False, pretty=True) is None

    model = load_json_model(json_path)
    # Because of changes to eq, to compare using the old format,
    # we need annotation.metadata.
    # TODO: get comments from cdiener
    # assert model.annotation.metadata == {
    #     "bigg.model": ["e_coli_core"],
    #     "doi": ["10.1128/ecosalplus.10.2.1"],
    #     "eco": ["ECO:0000004"],
    #     "ncbiprotein": ["16128336"],
    #     "pubmed": ["1111111"],
    #     "taxonomy": ["511145"],
    # }
    assert model.metadata.standardized == StandardizedAnnotationStore.from_data(
        ECOLI_MODEL_ANNOTATIONS
    )
    assert model.metadata.standardized == ECOLI_MODEL_ANNOTATIONS

    for met_id in model.metabolites.list_attr("id"):
        original_met_annot = annotation_model.metabolites.get_by_id(met_id).metadata
        new_met_annot = model.metabolites.get_by_id(met_id).metadata
        assert original_met_annot == new_met_annot

    for rxn_id in model.reactions.list_attr("id"):
        original_rxn_annot = annotation_model.reactions.get_by_id(rxn_id).metadata
        new_rxn_annot = model.reactions.get_by_id(rxn_id).metadata
        assert original_rxn_annot == new_rxn_annot


def test_read_write_sbml(annotation_model: Model, tmp_path: Path):
    """Test annotation consistency when writing and reading an SBML file."""
    out_path = tmp_path / "e_coli_core_json_writing.sbml"
    assert write_sbml_model(annotation_model, str(out_path)) is None

    model = read_sbml_model(str(out_path))
    # Because of changes to eq, to compare using the old format,
    # we need annotation.metadata
    # TODO: get comments from cdiener
    # assert model.annotation.metadata == {
    #     "bigg.model": ["e_coli_core"],
    #     "doi": ["10.1128/ecosalplus.10.2.1"],
    #     "eco": ["ECO:0000004"],
    #     "ncbiprotein": ["16128336"],
    #     "pubmed": ["1111111"],
    #     "taxonomy": ["511145"],
    # }
    assert model.metadata.standardized == StandardizedAnnotationStore.from_data(
        ECOLI_MODEL_ANNOTATIONS
    )
    assert model.metadata.standardized == ECOLI_MODEL_ANNOTATIONS

    for met_id in model.metabolites.list_attr("id"):
        original_met_annot = annotation_model.metabolites.get_by_id(met_id).metadata
        new_met_annot = model.metabolites.get_by_id(met_id).metadata
        assert original_met_annot == new_met_annot

    for rxn_id in model.reactions.list_attr("id"):
        original_rxn_annot = annotation_model.reactions.get_by_id(rxn_id).metadata
        new_rxn_annot = model.reactions.get_by_id(rxn_id).metadata
        assert original_rxn_annot == new_rxn_annot


def test_read_old_json_model(data_directory):
    """Test reading a schema v1 json model with old-style annotations."""
    model = load_json_model(Path(data_directory / "valid_annotation_format.json"))
    meta = model.metabolites[0]

    assert meta.annotation["bigg.reaction"] == ["PFK26"]
    assert meta.annotation["kegg.reaction"] == ["R02732"]
    assert meta.annotation == {
        "bigg.reaction": ["PFK26"],
        "kegg.reaction": ["R02732"],
        "rhea": ["15656"],
    }


# TODO: Fix this test. The mini.json model is currently not in the old format.
# def test_read_old_json_model(data_directory):
#     """Test reading the annotations of an old format JSON model."""
#     model = load_json_model(Path(data_directory) / "mini.json")
#     meta = model.metabolites[0]
#     # assert meta.annotation == {
#     #     "bigg.metabolite": ["13dpg"],
#     #     "biocyc": ["DPG"],
#     #     "chebi": [
#     #         "CHEBI:11881",
#     #         "CHEBI:16001",
#     #         "CHEBI:1658",
#     #         "CHEBI:20189",
#     #         "CHEBI:57604",
#     #     ],
#     #     "hmdb": ["HMDB01270"],
#     #     "kegg.compound": ["C00236"],
#     #     "pubchem.substance": ["3535"],
#     #     "reactome": ["REACT_29800"],
#     #     "seed.compound": ["cpd00203"],
#     #     "unipathway.compound": ["UPC00236"],
#     # }
#
#     # testing standardized
#     expected_cvterms = StandardizedAnnotationStore.from_data(
#         [{"qualifier": "bqb_is", "resources": RESOURCE_LIST}]
#     )
#     assert meta.metadata.standardized == expected_cvterms
#     assert meta.metadata.standardized == [
#         {"qualifier": "bqb_is", "identifeirs": RESOURCE_LIST}
#     ]
#


def test_cvtermlist_query():
    """Test the query functionality of StandardizedAnnotationStore."""
    resources = RESOURCE_LIST
    resources.extend(
        [
            "https://identifiers.org/uniprot/P69905",
            "https://identifiers.org/uniprot/P68871",
            "https://identifiers.org/kegg.compound/C00032",
            "https://identifiers.org/chebi/CHEBI:17627",
            "https://identifiers.org/chebi/CHEBI:43215",
            "https://identifiers.org/CHebi/CHEBI:11881",
        ]
    )
    cvtermlist = StandardizedAnnotationStore()
    for i, res in enumerate(resources):
        cvtermlist.extend(
            [
                StandardizedAnnotation(
                    qualifier=list(Qualifier._value2member_map_)[i], resources=res
                )
            ]
        )

    cvtermlist.append(
        StandardizedAnnotation(
            resources=ECO_EXAMPLE,
            annotations=[
                StandardizedAnnotation(
                    qualifier=Qualifier.Biological_isDescribedBy,
                    resources=PUBMED_EXAMPLE,
                )
            ],
            qualifier=list(Qualifier._value2member_map_)[19],
        )
    )
    print(cvtermlist)
    # assert isinstance(
    #     cvtermlist.query(search_function="bqm", attribute="qualifier"),
    #     StandardizedAnnotationStore,
    # )
    # The result type is now just a list, since a StandardizedAnnotationStore should be
    # used only for actual sets of annotations. This should maybe be a frozen variant of
    # the StandardizedAnnotationStore, but for now it is a list.
    assert isinstance(
        cvtermlist.query(search_function="bqm", attribute="qualifier"),
        StandardizedAnnotationList,
    )
    assert len(cvtermlist.query(search_function="bqm", attribute="qualifier")) == 6
    assert (
        len(cvtermlist.query(search_function="Modelling", attribute="qualifier")) == 6
    )
    assert (
        len(
            cvtermlist.query(search_function="bqb_isDescribedBy", attribute="qualifier")
        )
        == 1
    )
    assert (
        len(
            cvtermlist.query(
                search_function="Biological_isDescribedBy", attribute="qualifier"
            )
        )
        == 1
    )
    assert (
        len(cvtermlist.query(search_function=r"bqm_is\S+", attribute="qualifier")) == 3
    )
    assert (
        len(
            cvtermlist.query(
                search_function=lambda x: list(Qualifier._value2member_map_).index(
                    x.value
                )
                > 18,
                attribute="qualifier",
            )
        )
        == 1
    )

    # assert (
    #     len(
    #         cvtermlist.query(
    #             search_function=lambda x: x.nested_data,
    #             attribute="external_resources"
    #         )
    #     )
    #     == 1
    # )
    # assert (
    #     len(cvtermlist.query(search_function="chebi", attribute="external_resources"))
    #     == 7
    # )
    #
    # assert len(cvtermlist.query(search_function="chebi", attribute="resources")) == 7
    # assert (
    #     len(cvtermlist.query(search_function=r"[cC][hH]EBI", attribute="resources"))
    #     == 8
    # )
    # assert len(cvtermlist.query(search_function="pubmed", attribute="resources")) == 1
    #
    assert (
        len(cvtermlist.query(search_function=lambda x: x.qualifier.value == "bqm_is"))
        == 1
    )
    assert (
        len(
            cvtermlist.query(
                search_function=lambda x: x.qualifier.name == "Modelling_is"
            )
        )
        == 1
    )

    assert len(cvtermlist.query(search_function="chebi")) == 7
    assert len(cvtermlist.query(search_function=r"bqm_is\S+")) == 3
