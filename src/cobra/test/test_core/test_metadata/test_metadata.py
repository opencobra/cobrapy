"""
Tests for the metadata structures
"""

import json
import os
from os.path import join
from pathlib import Path

from cobra.core.metadata import CVTerms
from cobra.core.species import Species
from cobra.io import load_json_model, read_sbml_model, save_json_model, write_sbml_model


ecoli_model_annotation = {
    "bqb_hasTaxon": [{"resources": ["http://identifiers.org/taxonomy/511145"]}],
    "bqm_is": [
        {
            "nested_data": {
                "bqb_isDescribedBy": [
                    {"resources": ["https://identifiers.org/pubmed/1111111"]},
                    {"resources": ["https://identifiers.org/eco/ECO:0000004"]},
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

    # setting annotation via old annotation format
    s.annotation["chebi"] = ["CHEBI:43215", "CHEBI:11881"]

    # checking old (fixed) annotation format
    assert s.annotation == {"chebi": ["CHEBI:43215", "CHEBI:11881"]}

    # checking new cvterms
    cvt = CVTerms(
        {
            "bqb_is": [
                {
                    "resources": [
                        "https://identifiers.org/chebi/CHEBI:43215",
                        "https://identifiers.org/chebi/CHEBI:11881",
                    ]
                }
            ]
        }
    )
    assert s.annotation.cvterms == cvt

    # adding an SBO term
    s.annotation["sbo"] = ["SBO:0000123"]
    assert "chebi" in s.annotation
    assert "sbo" in s.annotation
    assert s.annotation == {
        "chebi": ["CHEBI:43215", "CHEBI:11881"],
        "sbo": ["SBO:0000123"],
    }


def test_nested_annotation(data_directory):
    # testing via cvterms
    with open(join(data_directory, "cvterms_nested.json"), "r") as f_cvterms:
        cvterms_data = json.load(f_cvterms)

    s = Species()
    s.annotation.add_cvterms(cvterms_data)
    assert s.annotation == {
        "uniprot": ["P69905", "P68871", "P69905"],
        "kegg.compound": ["C00032"],
        "chebi": ["CHEBI:17627"],
    }
    # check cvterms
    main_cvt = CVTerms(
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
                    "bqb_isDescribedBy": [
                        {"resources": ["https://identifiers.org/pubmed/1111111"]},
                        {"resources": ["https://identifiers.org/eco/000000"]},
                    ],
                },
            ]
        }
    )
    nested_cvt = CVTerms(
        {
            "bqb_isDescribedBy": [
                {"resources": ["https://identifiers.org/pubmed/1111111"]},
                {"resources": ["https://identifiers.org/eco/000000"]},
            ]
        }
    )
    assert s.annotation.cvterms == main_cvt
    nested_data = s.annotation.cvterms["bqb_hasPart"][1].nested_data
    assert nested_data == nested_cvt


def _read_ecoli_annotation_model(data_directory):
    test_xml = os.path.join(data_directory, "e_coli_core_for_annotation.xml")
    model = read_sbml_model(test_xml)
    return model


def test_cvterms_from_ecoli_xml(data_directory):
    model = _read_ecoli_annotation_model(data_directory)
    qualifier_list = ["bqb_hasTaxon", "bqm_is", "bqm_isDescribedBy"]
    nested_cvt = CVTerms(
        {
            "bqb_isDescribedBy": [
                {"resources": ["https://identifiers.org/pubmed/1111111"]},
                {"resources": ["https://identifiers.org/eco/ECO:0000004"]},
            ]
        }
    )
    ecoli_model_cvterm = CVTerms(ecoli_model_annotation)
    for key in qualifier_list:
        assert key in model.annotation.cvterms
    assert model.annotation.cvterms == ecoli_model_cvterm
    assert len(model.annotation.cvterms["bqm_isDescribedBy"]) == 2
    nested_data = model.annotation.cvterms["bqm_is"][0].nested_data
    assert nested_data == nested_cvt

    # check backwards compatibility
    assert model.annotation == {
        "taxonomy": ["511145"],
        "bigg.model": ["e_coli_core"],
        "doi": ["10.1128/ecosalplus.10.2.1"],
        "ncbigi": ["gi:16128336"],
    }


def test_writing_xml(data_directory, tmp_path):
    model = _read_ecoli_annotation_model(data_directory)
    assert write_sbml_model(model, tmp_path, "e_coli_core_writing.xml") is None


def test_read_write_json(data_directory, tmp_path):
    model = _read_ecoli_annotation_model(data_directory)
    json_path = join(tmp_path, "e_coli_core_json_writing.json")
    assert save_json_model(model, json_path, sort=False, pretty=True) is None

    model = load_json_model(json_path)
    assert model.annotation == {
        "taxonomy": ["511145"],
        "bigg.model": ["e_coli_core"],
        "doi": ["10.1128/ecosalplus.10.2.1"],
        "ncbigi": ["gi:16128336"],
    }
    assert model.annotation.cvterms == CVTerms(ecoli_model_annotation)


def test_read_old_json_model(data_directory):
    model = load_json_model(Path(data_directory) / "mini.json")
    meta = model.metabolites[0]
    assert meta.annotation == {
        "bigg.metabolite": ["13dpg"],
        "biocyc": ["DPG"],
        "chebi": [
            "CHEBI:16001",
            "CHEBI:1658",
            "CHEBI:20189",
            "CHEBI:57604",
            "CHEBI:11881",
        ],
        "hmdb": ["HMDB01270"],
        "kegg.compound": ["C00236"],
        "pubchem.substance": ["3535"],
        "reactome": ["REACT_29800"],
        "seed.compound": ["cpd00203"],
        "unipathway.compound": ["UPC00236"],
    }

    # testing cvterms
    expected_cvterms = CVTerms(
        {
            "bqb_is": [
                {
                    "resources": [
                        "https://identifiers.org/bigg.metabolite/13dpg",
                        "https://identifiers.org/biocyc/DPG",
                        "https://identifiers.org/chebi/CHEBI:16001",
                        "https://identifiers.org/chebi/CHEBI:1658",
                        "https://identifiers.org/chebi/CHEBI:20189",
                        "https://identifiers.org/chebi/CHEBI:57604",
                        "https://identifiers.org/chebi/CHEBI:11881",
                        "https://identifiers.org/hmdb/HMDB01270",
                        "https://identifiers.org/kegg.compound/C00236",
                        "https://identifiers.org/pubchem.substance/3535",
                        "https://identifiers.org/reactome/REACT_29800",
                        "https://identifiers.org/seed.compound/cpd00203",
                        "https://identifiers.org/unipathway.compound/UPC00236",
                    ]
                }
            ]
        }
    )
    assert meta.annotation.cvterms == expected_cvterms
