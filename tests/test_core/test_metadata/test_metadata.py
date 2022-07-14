"""Tests for the metadata structures."""

import json
import os
from os.path import join
from pathlib import Path

import pytest

from cobra.core.metadata import CVTerm, CVTermList, ExternalResources, Qualifier
from cobra.core.species import Species
from cobra.io import load_json_model, read_sbml_model, save_json_model, write_sbml_model


PUBMED_EXAMPLE = "https://identifiers.org/pubmed/1111111"
ECO_EXAMPLE = "https://identifiers.org/eco/ECO:0000004"
RESOURCE_LIST = [
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
    assert s.annotation.standardized == CVTermList()
    assert not s.annotation.keys()
    assert s.annotation.keyvaluepairs == {}
    assert s.annotation.history.creators == []
    assert s.annotation.history.modified_dates == []

    # setting annotation via old annotation format
    s.annotation["chebi"] = ["CHEBI:43215", "CHEBI:11881"]

    assert s.annotation.standardized.resources == {
        "https://identifiers.org/chebi/CHEBI:43215",
        "https://identifiers.org/chebi/CHEBI:11881",
    }

    # checking old (fixed) annotation format
    assert s.annotation == {"chebi": sorted(["CHEBI:43215", "CHEBI:11881"])}

    # checking new standardized
    cvt = CVTermList(
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

    assert s.annotation.standardized == cvt

    # adding an SBO term
    s.annotation["sbo"] = ["SBO:0000123"]
    assert "chebi" in s.annotation
    assert "sbo" in s.annotation
    assert s.annotation == {
        "chebi": sorted(["CHEBI:43215", "CHEBI:11881"]),
        "sbo": ["SBO:0000123"],
    }


def test_nested_annotation(data_directory):
    # testing via standardized
    with open(join(data_directory, "cvterms_nested.json"), "r") as f_cvterms:
        cvterms_data = json.load(f_cvterms)

    s = Species()
    s.annotation.add_cvterms(cvterms_data)
    assert s.annotation == {
        "chebi": ["CHEBI:17627"],
        "eco": ["000000"],
        "kegg.compound": ["C00032"],
        "pubmed": ["1111111"],
        "uniprot": ["P68871", "P69905"],
    }
    # check standardized
    main_cvt = {
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
    nested_cvt = {
        "bqb_isDescribedBy": [
            {"resources": [PUBMED_EXAMPLE]},
            {"resources": ["https://identifiers.org/eco/000000"]},
        ]
    }
    assert s.annotation.standardized == main_cvt
    nested_data = s.annotation.standardized[1].external_resources.nested_data
    assert nested_data == nested_cvt


def _read_ecoli_annotation_model(data_directory):
    test_xml = os.path.join(data_directory, "e_coli_core_for_annotation.xml")
    model = read_sbml_model(test_xml)
    return model


def test_cvterms_from_ecoli_xml(data_directory):
    model = _read_ecoli_annotation_model(data_directory)
    qualifier_set = {"bqb_hasTaxon", "bqm_is", "bqm_isDescribedBy"}
    nested_cvt = {
        "bqb_isDescribedBy": [
            {"resources": [PUBMED_EXAMPLE]},
            {"resources": [ECO_EXAMPLE]},
        ]
    }
    ecoli_model_cvterm = CVTermList.from_dict(ecoli_model_annotation)
    xml_model_cvterms = model.annotation.standardized
    model_cvterms_qualifier_set = {qual.name for qual in xml_model_cvterms.qualifiers}
    assert qualifier_set == model_cvterms_qualifier_set
    assert xml_model_cvterms == ecoli_model_cvterm
    assert (
        len(model.annotation.standardized.query("bqm_isDescribedBy", "qualifier")) == 2
    )
    nested_data = model.annotation.standardized.query("bqm_is", "qualifier")[
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
    assert (
        write_sbml_model(model, str(tmp_path.joinpath("e_coli_core_writing.xml")))
        is None
    )


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
    assert model.annotation.standardized == CVTermList.from_dict(ecoli_model_annotation)
    assert model.annotation.standardized == ecoli_model_annotation


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

    # testing standardized
    expected_cvterms = CVTermList.from_dict({"bqb_is": [{"resources": RESOURCE_LIST}]})
    assert meta.annotation.standardized == expected_cvterms
    assert meta.annotation.standardized == {"bqb_is": [{"resources": RESOURCE_LIST}]}


def test_cvtermlist_query():
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
    cvtermlist = CVTermList()
    for i, res in enumerate(resources):
        cvtermlist.extend([CVTerm(qualifier=Qualifier(i), ex_res=resources[i])])

    cvtermlist.append(
        CVTerm(
            ex_res={
                "resources": ECO_EXAMPLE,
                "nested_data": CVTerm(
                    qualifier="bqb_isDescribedBy", ex_res=PUBMED_EXAMPLE
                ),
            },
            qualifier=Qualifier(19),
        )
    )
    assert isinstance(
        cvtermlist.query(search_function="bqm", attribute="qualifier"), CVTermList
    )
    assert len(cvtermlist.query(search_function="bqm", attribute="qualifier")) == 6
    assert (
        len(cvtermlist.query(search_function=r"bqm_is\S+", attribute="qualifier")) == 3
    )
    assert (
        len(
            cvtermlist.query(
                search_function=lambda x: x.value > 18, attribute="qualifier"
            )
        )
        == 1
    )

    assert (
        len(
            cvtermlist.query(
                search_function=lambda x: x.nested_data, attribute="external_resources"
            )
        )
        == 1
    )
    with pytest.raises(TypeError):
        assert (
            len(
                cvtermlist.query(
                    search_function="chebi", attribute="external_resources"
                )
            )
            == 1
        )

    assert len(cvtermlist.query(search_function="chebi", attribute="resources")) == 7
    assert (
        len(cvtermlist.query(search_function=r"[cC][hH]EBI", attribute="resources"))
        == 8
    )
    assert len(cvtermlist.query(search_function="pubmed", attribute="resources")) == 1
