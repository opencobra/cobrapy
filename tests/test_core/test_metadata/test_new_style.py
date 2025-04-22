"""Tests for the metadata structures."""

import json
from pathlib import Path
from pprint import pprint

import pytest

from cobra import Model
from cobra.core.metadata import Qualifier, StandardizedAnnotation
from cobra.core.metadata.custom import CustomAnnotation
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
        "resources": ["http://identifiers.org/tanomy/511145"],
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


def test_annotation() -> None:
    """Test creating an annotation manually."""
    s = Species()

    # assert s.annotations is None
    annotation_1 = StandardizedAnnotation("https://identifiers.org/go/GO:0007268")
    # Default qualifier should be bqb_is
    assert annotation_1.qualifier == Qualifier.Biological_is
    annotation_2 = StandardizedAnnotation(
        [
            "https://identifiers.org/wikipathways/WP179",
            "https://identifiers.org/reactome/REACT_152",
        ],
        qualifier=Qualifier.Biological_isVersionOf,
    )
    assert len(annotation_2.resources) == 2
    nested_annotations = [
        StandardizedAnnotation(
            "https://identifiers.org/pubmed/1111111",
            qualifier=Qualifier.Biological_isDescribedBy,
        ),
        StandardizedAnnotation(
            "https://identifiers.org/eco/ECO:0000004",
            qualifier=Qualifier.Biological_isDescribedBy,
        ),
    ]

    annotation_3 = StandardizedAnnotation(
        qualifier=Qualifier.Biological_occursIn,
        resources=["https://identifiers.org/go/GO:0005764"],
        annotations=nested_annotations,
    )
    assert len(annotation_3.annotations) == 2
    s.add_annotations([annotation_1, annotation_2, annotation_3])
    # assert len(s.annotations) == 3
    assert len(s.metadata.standardized) == 3
    s.remove_annotations([annotation_2])
    assert len(s.metadata.standardized) == 2
    s.add_annotations(annotation_2)
    assert len(s.metadata.standardized) == 3
    annotation_2.remove_from_parent()
    assert len(s.metadata.standardized) == 2
    # assert s.annotations.custom is None
    # assert (
    #     len(
    #         s.annotations.standardized.by_qualifier(Qualifier.Biological_is)[
    #             0
    #         ].resources
    #     )
    #     == 2
    # )
    custom_1 = CustomAnnotation(key="cobra_flag", value="starred")
    custom_2 = CustomAnnotation(key="cobra_url", uri=COBRA_URL)
    s.add_annotations([custom_1, custom_2])
    assert len(s.metadata.custom) == 2
    s.remove_annotations(custom_1)
    assert len(s.metadata.custom) == 1
    assert s.metadata.custom["cobra_url"].uri == COBRA_URL


def test_read_write_sbml(annotation_model: Model, tmp_path: Path):
    """Test annotation consistency when writing and reading an SBML file."""
    out_path = tmp_path / "e_coli_core_json_writing.sbml"
    assert write_sbml_model(annotation_model, str(out_path)) is None

    model = read_sbml_model(str(out_path))
    print(model.metadata)
    print(model.metadata.standardized)
    for ann in model.metadata.standardized:
        print(ann)
    print(model.metadata.standardized.resources)
    pprint(model.metadata.standardized.to_records())
    pprint(model.metadata.standardized.to_list_of_dicts())
    assert len(model.metadata.standardized) == 4
    ann_dict = [
        {
            "resources": [
                "http://identifiers.org/taxonomy/511145",
            ],
            "qualifier": "bqb_hasTaxon",
        },
        {
            "annotations": [
                {
                    "resources": [
                        "https://identifiers.org/pubmed/1111111",
                    ],
                    "qualifier": "bqb_isDescribedBy",
                },
                {
                    "resources": [
                        "https://identifiers.org/eco/ECO:0000004",
                    ],
                    "qualifier": "bqb_isDescribedBy",
                },
            ],
            "resources": [
                "http://identifiers.org/bigg.model/e_coli_core",
            ],
            "qualifier": "bqm_is",
        },
        {
            "resources": [
                "http://identifiers.org/doi/10.1128/ecosalplus.10.2.1",
            ],
            "qualifier": "bqm_isDescribedBy",
        },
        {
            "resources": [
                "http://identifiers.org/ncbiprotein/16128336",
            ],
            "qualifier": "bqm_isDescribedBy",
        },
    ]
    assert model.metadata.standardized == ann_dict
    # Because of changes to eq, to compare using the old format,
    # we need annotation.annotations
    # TODO: get comments from cdiener
    # assert model.annotation.annotations == {
    #     "bigg.model": ["e_coli_core"],
    #     "doi": ["10.1128/ecosalplus.10.2.1"],
    #     "eco": ["ECO:0000004"],
    #     "ncbiprotein": ["16128336"],
    #     "pubmed": ["1111111"],
    #     "taxonomy": ["511145"],
    # }
    # assert model.annotation.standardized == CVTermList.from_data(
    #     ECOLI_MODEL_ANNOTATIONS
    # )
    # assert model.annotation.standardized == ECOLI_MODEL_ANNOTATIONS
    #
    # for met_id in model.metabolites.list_attr("id"):
    #     original_met_annot = annotation_model.metabolites.get_by_id(met_id).annotation
    #     new_met_annot = model.metabolites.get_by_id(met_id).annotation
    #     assert original_met_annot == new_met_annot
    #
    # for rxn_id in model.reactions.list_attr("id"):
    #     original_rxn_annot = annotation_model.reactions.get_by_id(rxn_id).annotation
    #     new_rxn_annot = model.reactions.get_by_id(rxn_id).annotation
    #     assert original_rxn_annot == new_rxn_annot
