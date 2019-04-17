from os.path import join

from cobra.io import read_sbml_model, write_sbml_model


def _check_sbml_annotations(model):
    """Checks the annotations from the annotation.xml. """
    assert model is not None

    # model annotation
    # {'bigg.model': 'e_coli_core', 'doi': '10.1128/ecosalplus.10.2.1',
    #  'taxonomy': '511145'}
    annotation = model.annotation
    assert annotation is not None
    assert len(annotation) == 3
    for key in ["bigg.model", "doi", "taxonomy"]:
        assert key in annotation
    assert annotation["bigg.model"] == "e_coli_core"
    assert annotation["doi"] == "10.1128/ecosalplus.10.2.1"
    assert annotation["taxonomy"] == "511145"

    # gene annotation
    # {'asap': 'ABE-0006162', 'ncbigene': '946368', 'uniprot': 'P33221',
    #  'ncbigi': 'gi:16129802', 'ecogene': 'EG11809'}
    annotation = model.genes.G1.annotation
    assert len(annotation) == 5
    for key in ["asap", "ncbigene", "uniprot", "ncbigi", "ecogene"]:
        assert key in annotation
    assert annotation["asap"] == "ABE-0006162"
    assert annotation["ncbigene"] == "946368"
    assert annotation["uniprot"] == "P33221"
    assert annotation["ncbigi"] == "gi:16129802"
    assert annotation["ecogene"] == "EG11809"

    # compartment annotation
    # FIXME: add tests with first class compartment model
    # annotation = model.compartments.c.annotation
    # assert len(annotation) == 1
    # for key in ["bigg.compartment"]:
    #    assert key in annotation
    # assert annotation["bigg.compartment"] == "c"

    # metabolite/species annotation
    # {'inchi': 'InChI=1S/C3H8O2/c1-3(5)2-4/h3-5H,2H2,1H3/t3-/m0/s1',
    # 'bigg.metabolite': '13dpg', 'chebi': ['CHEBI:11881', 'CHEBI:16001',
    # 'CHEBI:1658', 'CHEBI:20189', 'CHEBI:57604'],
    # 'metanetx.chemical': 'MNXM261',
    # 'kegg.compound': ['C00236', 'C02917'],
    # 'seed.compound': 'cpd00203', 'hmdb': ['HMDB62758', 'HMDB06213'],
    # 'biocyc': 'META:DPG'}
    annotation = model.metabolites.A.annotation
    for key in ["inchi", "bigg.metabolite", "chebi", "metanetx.chemical",
                "kegg.compound", "seed.compound", "hmdb", "biocyc"]:
        assert key in annotation
    assert annotation[
               "inchi"] == "InChI=1S/C3H8O2/c1-3(5)2-4/h3-5H,2H2,1H3/t3-/m0/s1"  # noqa: E501

    # reaction annotation
    # {'kegg.reaction': 'R00228', 'sbo': 'SBO:0000375',
    # 'ec-code': '1.2.1.10', 'rhea': ['23288', '23289', '23290', '23291'],
    # 'metanetx.reaction': 'MNXR95210', 'bigg.reaction': 'ACALD',
    # 'biocyc': 'META:ACETALD-DEHYDROG-RXN'}
    annotation = model.reactions.R1.annotation
    for key in ["kegg.reaction", "sbo", "ec-code", "rhea",
                "metanetx.reaction", "bigg.reaction", "biocyc"]:
        assert key in annotation
    assert annotation["biocyc"] == 'META:ACETALD-DEHYDROG-RXN'


def test_read_sbml_annotations(data_directory):
    """Test reading and writing annotations."""
    with open(join(data_directory, "annotation.xml"), "r") as f_in:
        model1 = read_sbml_model(f_in)
        _check_sbml_annotations(model1)


def test_read_write_sbml_annotations(data_directory, tmp_path):
    """Test reading and writing annotations."""
    with open(join(data_directory, "annotation.xml"), "r") as f_in:
        model1 = read_sbml_model(f_in)

    sbml_path = join(str(tmp_path), "test.xml")
    with open(sbml_path, "w") as f_out:
        write_sbml_model(model1, f_out)

    with open(sbml_path, "r") as f_in:
        model2 = read_sbml_model(f_in)
        _check_sbml_annotations(model2)
