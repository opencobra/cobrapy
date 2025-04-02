"""Tests for the metadata structures."""

import json
from pathlib import Path

from cobra.core.metadata.keyvaluepairs import CustomAnnotation
import pytest

from cobra.core.metadata import StandardizedAnnotation, Qualifier
from cobra.core.species import Species

# from cobra.io import load_json_model, read_sbml_model, save_json_model, write_sbml_model


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
        "external_resources": {"resources": ["http://identifiers.org/taxonomy/511145"]},
    },
    {
        "qualifier": "bqm_is",
        "external_resources": {
            "resources": ["http://identifiers.org/bigg.model/e_coli_core"],
            "nested_data": [
                {
                    "qualifier": "bqb_isDescribedBy",
                    "external_resources": {"resources": [PUBMED_EXAMPLE]},
                },
                {
                    "qualifier": "bqb_isDescribedBy",
                    "external_resources": {"resources": [ECO_EXAMPLE]},
                },
            ],
        },
    },
    {
        "qualifier": "bqm_isDescribedBy",
        "external_resources": {
            "resources": ["http://identifiers.org/doi/10.1128/ecosalplus.10.2.1"]
        },
    },
    {
        "qualifier": "bqm_isDescribedBy",
        "external_resources": {
            "resources": ["http://identifiers.org/ncbiprotein/16128336"]
        },
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
    assert len(annotation_2.identifiers) == 2
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
        identifiers=["https://identifiers.org/go/GO:0005764"],
        annotations=nested_annotations,
    )
    assert len(annotation_3.annotations) == 2
    s.add_annotations([annotation_1, annotation_2, annotation_3])
    # assert len(s.annotations) == 3
    assert len(s.annotations.standardized) == 3
    s.remove_annotations([annotation_2])
    assert len(s.annotations.standardized) == 2
    s.add_annotations(annotation_2)
    assert len(s.annotations.standardized) == 3
    annotation_2.remove_from_object()
    assert len(s.annotations.standardized) == 2
    # assert s.annotations.custom is None
    # assert (
    #     len(
    #         s.annotations.standardized.by_qualifier(Qualifier.Biological_is)[
    #             0
    #         ].identifiers
    #     )
    #     == 2
    # )
    custom_1 = CustomAnnotation(key="cobra_flag", value="starred")
    custom_2 = CustomAnnotation(key="cobra_url", uri=COBRA_URL)
    s.add_annotations([custom_1, custom_2])
    assert len(s.annotations.custom) == 2
    s.remove_annotations(custom_1)
    assert len(s.annotations.custom) == 1
    assert s.annotations.custom["cobra_url"].uri == COBRA_URL
