from cobra.core.species import Species
from cobra.core.metadata import *
import pytest
import json
from pathlib import Path
from cobra.io import read_sbml_model, write_sbml_model, load_json_model, save_json_model

def test_annotation_working():
    # a cobra component
    s = Species()
    assert s.annotation == {} # nothing set for annotation, so empty dict

    # setting annotation via old format
    s.annotation["chebi"] = ["CHEBI:43215", "CHEBI:11881"]
    assert s.annotation == {"chebi": ["CHEBI:43215", "CHEBI:11881"]}
    cvt = CVTerms({'bqb_is': [
                        {'resources': ['https://identifiers.org/chebi/CHEBI:43215',
                                       'https://identifiers.org/chebi/CHEBI:11881']
                        }
                    ]
                })
    assert str(s.annotation.cvterms) == str(cvt) == "{'bqb_is': [{'resources': ['https://identifiers.org/chebi/CHEBI:43215', 'https://identifiers.org/chebi/CHEBI:11881']}]}"

    # testing sbo term
    s.annotation["sbo"] = "SBO:0000123"
    assert s.annotation == {'chebi':
                                ['CHEBI:43215',
                                 'CHEBI:11881'],
                            'sbo': 'SBO:0000123'
                           }

    assert "chebi" in s.annotation
    assert "sbo" in s.annotation

    # testing via cvterms
    with open(Path(__file__).parent / f"cvterms_nested.json", "r") as f_cvterms:
        cvterms_data = json.load(f_cvterms)
    s.annotation.add_cvterms(cvterms_data)
    assert s.annotation == {'chebi':
                                ['CHEBI:43215',
                                 'CHEBI:11881',
                                 'CHEBI:17627'],
                            'sbo': 'SBO:0000123',
                            'uniprot':
                                ['P69905',
                                 'P68871',
                                 'P69905'],
                            'kegg.compound':
                                ['C00032']
                            }
    assert str(s.annotation.cvterms) == "{'bqb_is': [{'resources': ['https://identifiers.org/chebi/CHEBI:43215', 'https://identifiers.org/chebi/CHEBI:11881']}], 'bqb_hasPart': [{'resources': ['https://identifiers.org/uniprot/P69905', 'https://identifiers.org/uniprot/P68871', 'https://identifiers.org/kegg.compound/C00032']}, {'resources': ['https://identifiers.org/uniprot/P69905', 'https://www.uniprot.org/uniprot/P68871', 'https://identifiers.org/chebi/CHEBI:17627'], 'nested_data': {'bqb_isDescribedBy': [{'resources': ['https://identifiers.org/pubmed/1111111']}, {'resources': ['https://identifiers.org/eco/000000']}]}}]}"


def test_reading_xml():
    model = read_sbml_model("cobra/test/data/e_coli_core_for_annotation.xml")
    assert model.annotation == {
            'taxonomy': ['511145'],
            'bigg.model': ['e_coli_core'],
            'doi': ['10.1128/ecosalplus.10.2.1'],
            'ncbigi': ['gi:16128336']
        }
    assert str(model.annotation.cvterms) == "{'bqb_hasTaxon': [{'resources': ['http://identifiers.org/taxonomy/511145']}], 'bqm_is': [{'resources': ['http://identifiers.org/bigg.model/e_coli_core'], 'nested_data': {'bqb_isDescribedBy': [{'resources': ['https://identifiers.org/pubmed/1111111']}, {'resources': ['https://identifiers.org/eco/ECO:0000004']}]}}], 'bqm_isDescribedBy': [{'resources': ['http://identifiers.org/doi/10.1128/ecosalplus.10.2.1']}, {'resources': ['http://identifiers.org/ncbigi/gi:16128336']}]}"
    assert str(model.annotation.history) == "{'creators': [{'first_name': 'Matthias', 'last_name': 'Koenig', 'email': 'koenigmx@hu-berlin.de', 'organization_name': 'Humboldt-University Berlin, Institute for Theoretical Biology'}], 'created': '2019-03-06T14:40:55Z', 'modified': ['2019-03-06T14:40:55Z']}"


def test_writing_xml():
    model = read_sbml_model("cobra/test/data/e_coli_core_for_annotation.xml")
    assert write_sbml_model(model, "cobra/core/metadata/examples/e_coli_core_writing.xml") is None


def test_write_json():
    model = read_sbml_model("cobra/test/data/e_coli_core_for_annotation.xml")
    assert save_json_model(model, "cobra/core/metadata/examples/e_coli_core_json_writing.json", False, True) is None


def test_read_new_json_model():
    model = load_json_model("cobra/core/metadata/examples/e_coli_core_json_writing.json")
    assert model.annotation == {
            'taxonomy': ['511145'],
            'bigg.model': ['e_coli_core'],
            'doi': ['10.1128/ecosalplus.10.2.1'],
            'ncbigi': ['gi:16128336']
        }
    assert str(model.annotation.cvterms) == "{'bqb_hasTaxon': [{'resources': ['http://identifiers.org/taxonomy/511145']}], 'bqm_is': [{'resources': ['http://identifiers.org/bigg.model/e_coli_core'], 'nested_data': {'bqb_isDescribedBy': [{'resources': ['https://identifiers.org/pubmed/1111111']}, {'resources': ['https://identifiers.org/eco/ECO:0000004']}]}}], 'bqm_isDescribedBy': [{'resources': ['http://identifiers.org/doi/10.1128/ecosalplus.10.2.1']}, {'resources': ['http://identifiers.org/ncbigi/gi:16128336']}]}"
    assert str(model.annotation.history) == "{'creators': [{'first_name': 'Matthias', 'last_name': 'Koenig', 'email': 'koenigmx@hu-berlin.de', 'organization_name': 'Humboldt-University Berlin, Institute for Theoretical Biology'}], 'created': '2019-03-06T14:40:55Z', 'modified': ['2019-03-06T14:40:55Z']}"


def test_read_old_json_model():
    model = load_json_model("cobra/test/data/mini.json")
    meta = model.metabolites[0]
    assert meta.annotation == {
        'bigg.metabolite': ['13dpg'],
        'biocyc': ['DPG'],
        'chebi': ['CHEBI:16001',
                  'CHEBI:1658',
                  'CHEBI:20189',
                  'CHEBI:57604',
                  'CHEBI:11881'],
        'hmdb': ['HMDB01270'],
        'kegg.compound': ['C00236'],
        'pubchem.substance': ['3535'],
        'reactome': ['REACT_29800'],
        'seed.compound': ['cpd00203'],
        'unipathway.compound': ['UPC00236']
    }
    assert str(meta.annotation.cvterms) == "{'bqb_is': [{'resources': ['https://identifiers.org/bigg.metabolite/13dpg', 'https://identifiers.org/biocyc/DPG', 'https://identifiers.org/chebi/CHEBI:16001', 'https://identifiers.org/chebi/CHEBI:1658', 'https://identifiers.org/chebi/CHEBI:20189', 'https://identifiers.org/chebi/CHEBI:57604', 'https://identifiers.org/chebi/CHEBI:11881', 'https://identifiers.org/hmdb/HMDB01270', 'https://identifiers.org/kegg.compound/C00236', 'https://identifiers.org/pubchem.substance/3535', 'https://identifiers.org/reactome/REACT_29800', 'https://identifiers.org/seed.compound/cpd00203', 'https://identifiers.org/unipathway.compound/UPC00236']}]}"
