"""Tests for the metadata structures."""

import json
import os
from os.path import join
from pathlib import Path

from cobra.core.metadata import CVTerm, CVTerms, ExternalResources
from cobra.core.species import Species
from cobra.io import load_json_model, read_sbml_model, save_json_model, write_sbml_model

PUBMED_EXAMPLE = "https://identifiers.org/pubmed/1111111"
ECO_EXAMPLE = "https://identifiers.org/eco/ECO:0000004"

ecoli_model_annotation = {
    "bqb_hasTaxon": [{"resources": ["http://identifiers.org/taxonomy/511145"]}],
    "bqm_is": [
        {
            "nested_data": {
                "bqb_isDescribedBy": [
                    {"resources": [PUBMED_EXAMPLE]},
                    {"resources": [ECO_EXAMPLE]},
                ]
            },
            "resources": ["http://identifiers.org/bigg.model/e_coli_core"],
        }
    ],
    "bqm_isDescribedBy": [
        {"resources": ["http://identifiers.org/doi/10.1128/ecosalplus.10.2.1"]},
        {"resources": ["http://identifiers.org/ncbigi/gi:16128336"]},
    ],
}


def test_annotation():
    # a cobra component
    s = Species()
    assert s.annotation == {}  # nothing set for annotation, so empty dict
    assert s.annotation.cvterms == CVTerms()
    assert not s.annotation.keys()
    assert s.annotation.keyvaluepairs == {}
    assert s.annotation.history.creators == []
    assert s.annotation.history.modified_dates == []

    # setting annotation via old annotation format
    s.annotation["chebi"] = ["CHEBI:43215", "CHEBI:11881"]

    assert s.annotation.cvterms.resources == {
        "https://identifiers.org/chebi/CHEBI:43215",
        "https://identifiers.org/chebi/CHEBI:11881",
    }

    # checking old (fixed) annotation format
    assert s.annotation == {"chebi": sorted(["CHEBI:43215", "CHEBI:11881"])}

    # checking new cvterms
    cvt = CVTerms(
        [
            CVTerm(
                qualifier="bqb_is",
                ex_res=ExternalResources(
                    resources=["https://identifiers.org/chebi/CHEBI:43215"]
                ),
            ),
            CVTerm(
                qualifier="bqb_is",
                ex_res=ExternalResources(
                    resources=["https://identifiers.org/chebi/CHEBI:11881"]
                ),
            ),
        ]
    )

    assert s.annotation.cvterms == cvt

    # adding an SBO term
    s.annotation["sbo"] = ["SBO:0000123"]
    assert "chebi" in s.annotation
    assert "sbo" in s.annotation
    assert s.annotation == {
        "chebi": sorted(["CHEBI:43215", "CHEBI:11881"]),
        "sbo": ["SBO:0000123"],
    }


def test_nested_annotation(data_directory):
    # testing via cvterms
    with open(join(data_directory, "cvterms_nested.json"), "r") as f_cvterms:
        cvterms_data = json.load(f_cvterms)

    s = Species()
    s.annotation.add_cvterms(CVTerms.from_dict(cvterms_data))
    assert s.annotation == {
        "chebi": ["CHEBI:17627"],
        "eco": ["000000"],
        "kegg.compound": ["C00032"],
        "pubmed": ["1111111"],
        "uniprot": ["P68871", "P69905"],
    }
    # check cvterms
    main_cvt = CVTerms.from_dict(
        {
            "bqb_hasPart": [
                {
                    "resources": [
                        "https://identifiers.org/uniprot/P69905",
                        "https://identifiers.org/uniprot/P68871",
                        "https://identifiers.org/kegg.compound/C00032",
                    ]
                },
                {
                    "resources": [
                        "https://identifiers.org/uniprot/P69905",
                        "https://www.uniprot.org/uniprot/P68871",
                        "https://identifiers.org/chebi/CHEBI:17627",
                    ],
                    "nested_data": {
                        "bqb_isDescribedBy": [
                            {"resources": [PUBMED_EXAMPLE]},
                            {"resources": ["https://identifiers.org/eco/000000"]},
                        ],
                    },
                },
            ]
        }
    )
    nested_cvt = CVTerms.from_dict(
        {
            "bqb_isDescribedBy": [
                {"resources": [PUBMED_EXAMPLE]},
                {"resources": ["https://identifiers.org/eco/000000"]},
            ]
        }
    )
    assert s.annotation.cvterms == main_cvt
    nested_data = s.annotation.cvterms.query_qualifier("bqb_hasPart")[
        1
    ].external_resources.nested_data
    assert nested_data == nested_cvt


def _read_ecoli_annotation_model(data_directory):
    test_xml = os.path.join(data_directory, "e_coli_core_for_annotation.xml")
    model = read_sbml_model(test_xml)
    return model


def test_cvterms_from_ecoli_xml(data_directory):
    model = _read_ecoli_annotation_model(data_directory)
    qualifier_set = {"bqb_hasTaxon", "bqm_is", "bqm_isDescribedBy"}
    nested_cvt = CVTerms.from_dict(
        {
            "bqb_isDescribedBy": [
                {"resources": [PUBMED_EXAMPLE]},
                {"resources": [ECO_EXAMPLE]},
            ]
        }
    )
    ecoli_model_cvterm = CVTerms.from_dict(ecoli_model_annotation)
    xml_model_cvterms = model.annotation.cvterms
    model_cvterms_qualifier_set = {qual.name for qual in xml_model_cvterms.qualifiers}
    assert qualifier_set == model_cvterms_qualifier_set
    assert xml_model_cvterms == ecoli_model_cvterm
    assert len(xml_model_cvterms.query_qualifier("bqm_isDescribedBy")) == 2
    nested_data = xml_model_cvterms.query_qualifier("bqm_is")[
        0
    ].external_resources.nested_data
    assert nested_data == nested_cvt

    # check backwards compatibility
    assert model.annotation == {
        "bigg.model": ["e_coli_core"],
        "doi": ["10.1128/ecosalplus.10.2.1"],
        "eco": ["ECO:0000004"],
        "ncbigi": ["gi:16128336"],
        "pubmed": ["1111111"],
        "taxonomy": ["511145"],
    }


def test_writing_xml(data_directory, tmp_path):
    model = _read_ecoli_annotation_model(data_directory)
    assert write_sbml_model(model, tmp_path.joinpath("e_coli_core_writing.xml")) is None


def test_read_write_json(data_directory, tmp_path):
    model = _read_ecoli_annotation_model(data_directory)
    json_path = join(tmp_path, "e_coli_core_json_writing.json")
    assert save_json_model(model, json_path, sort=False, pretty=True) is None

    model = load_json_model(json_path)
    assert model.annotation == {
        "bigg.model": ["e_coli_core"],
        "doi": ["10.1128/ecosalplus.10.2.1"],
        "eco": ["ECO:0000004"],
        "ncbigi": ["gi:16128336"],
        "pubmed": ["1111111"],
        "taxonomy": ["511145"],
    }
    assert model.annotation.cvterms == CVTerms.from_dict(ecoli_model_annotation)


def test_read_old_json_model(data_directory):
    model = load_json_model(Path(data_directory) / "mini.json")
    meta = model.metabolites[0]
    assert meta.annotation == {
        "bigg.metabolite": ["13dpg"],
        "biocyc": ["DPG"],
        "chebi": [
            "CHEBI:11881",
            "CHEBI:16001",
            "CHEBI:1658",
            "CHEBI:20189",
            "CHEBI:57604",
        ],
        "hmdb": ["HMDB01270"],
        "kegg.compound": ["C00236"],
        "pubchem.substance": ["3535"],
        "reactome": ["REACT_29800"],
        "seed.compound": ["cpd00203"],
        "unipathway.compound": ["UPC00236"],
    }

    # testing cvterms
    expected_cvterms = CVTerms.from_dict(
        {
            "bqb_is": [
                {
                    "resources": [
                        "http://identifiers.org/bigg.metabolite/13dpg",
                        "http://identifiers.org/biocyc/DPG",
                        "http://identifiers.org/chebi/CHEBI:16001",
                        "http://identifiers.org/chebi/CHEBI:1658",
                        "http://identifiers.org/chebi/CHEBI:20189",
                        "http://identifiers.org/chebi/CHEBI:57604",
                        "http://identifiers.org/chebi/CHEBI:11881",
                        "http://identifiers.org/hmdb/HMDB01270",
                        "http://identifiers.org/kegg.compound/C00236",
                        "http://identifiers.org/pubchem.substance/3535",
                        "http://identifiers.org/reactome/REACT_29800",
                        "http://identifiers.org/seed.compound/cpd00203",
                        "http://identifiers.org/unipathway.compound/UPC00236",
                    ]
                }
            ]
        }
    )
    assert meta.annotation.cvterms == expected_cvterms
