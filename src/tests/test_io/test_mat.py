"""Test functionalities of I/O in MATLAB (.mat) format."""

from os.path import join
from pickle import load
from typing import TYPE_CHECKING

import pytest

from cobra import io


try:
    import scipy
except ImportError:
    scipy = None


if TYPE_CHECKING:
    import py.path

    from cobra import Model
    from typing import Callable


@pytest.fixture(scope="function")
def raven_model(data_directory: str) -> "Model":
    """Fixture for RAVEN model."""
    with open(join(data_directory, "raven.pickle"), "rb") as infile:
        return load(infile)


@pytest.mark.skipif(scipy is None, reason="scipy unavailable")
# @pytest.mark.parametrize("ref_model, filename",
#                          [(pytest.fixture_request("mini_model"),
#                            "mini.mat"),
#                           (pytest.fixture_request("raven_model"),
#                            "raven.mat")])
# TODO: wait for pytest.fixture_request() to get approved
def test_load_matlab_model(
    compare_models: Callable,
    data_directory: str,
    mini_model: "Model",
    raven_model: "Model",
) -> None:
    """Test the reading of MAT model."""
    mini_mat_model = io.load_matlab_model(join(data_directory, "mini.mat"))
    raven_mat_model = io.load_matlab_model(join(data_directory, "raven.mat"))
    assert compare_models(mini_model, mini_mat_model) is None
    assert compare_models(raven_model, raven_mat_model) is None


# @pytest.mark.xfail(reason="localPath not supported yet")
@pytest.mark.skipif(scipy is None, reason="scipy unavailable")
# @pytest.mark.parametrize("model, filename",
#                          [(pytest.fixture_request("mini_model"),
#                            "mini.mat"),
#                           (pytest.fixture_request("raven_model"),
#                            "raven.mat")])
# TODO: wait for pytest.fixture_request() to get approved
def test_save_matlab_model(
    tmpdir: "py.path.local", mini_model: "Model", raven_model: "Model"
) -> None:
    """Test the writing of MAT model."""
    mini_output_file = tmpdir.join("mini.mat")
    raven_output_file = tmpdir.join("raven.mat")
    # scipy.io.savemat() doesn't support anything other than
    # str or file-stream object, hence the str conversion
    io.save_matlab_model(mini_model, str(mini_output_file))
    io.save_matlab_model(raven_model, str(raven_output_file))
    assert mini_output_file.check()
    assert raven_output_file.check()


@pytest.mark.skipif(scipy is None, reason="scipy unavailable")
def test_large_bounds(tmpdir: "py.path.local", model: "Model") -> None:
    """Verify that mat bounds don't get broken by the config defaults."""
    model.reactions[0].bounds = -1e6, 1e6
    filepath = str(tmpdir.join("model.mat"))
    io.save_matlab_model(model, filepath)
    read = io.load_matlab_model(filepath)
    assert read.reactions[0].bounds == (-1e6, 1e6)


@pytest.mark.skipif(scipy is None, reason="scipy unavailable")
def test_read_rewrite_matlab_model(
    tmpdir: "py.path.local", data_directory: str
) -> None:
    """Verify that rewritten matlab model is identical to original."""
    mini_mat_model = io.load_matlab_model(join(data_directory, "mini.mat"))
    raven_mat_model = io.load_matlab_model(join(data_directory, "raven.mat"))
    mini_output_file = tmpdir.join("mini.mat")
    raven_output_file = tmpdir.join("raven.mat")
    # scipy.io.savemat() doesn't support anything other than
    # str or file-stream object, hence the str conversion
    io.save_matlab_model(mini_mat_model, str(mini_output_file))
    io.save_matlab_model(raven_mat_model, str(raven_output_file))
    mini_mat_model_reload = io.load_matlab_model(str(mini_output_file))
    raven_mat_model_reload = io.load_matlab_model(str(raven_output_file))
    assert compare_models(mini_mat_model, mini_mat_model_reload) is None
    assert compare_models(raven_mat_model, raven_mat_model_reload) is None


def _fix_xml_annotation_to_identifiers(model: "Model") -> None:
    """Some XML models with cobra have annotations that do not match identifiers.org.

    This function will fix the dict keys of annotations to match identifiers.org.
    Eventually, the XML models should be fixed and cobrapy should be strict, but this is
    part of SBML rewriting of annotations
    see: https://github.com/opencobra/cobrapy/issues/684

    It also changes met formulas from empty string to None (which is the default
    when creating a metabolite with no fomula given) and strips spaces from reaction
    names.

    Parameters
    ----------
    model: Model
        A model to fix
    """
    for met in model.metabolites:
        if met.formula == "":
            met.formula = None
        if len(met.annotation):
            if "chebi" in met.annotation.keys():
                met.annotation["CHEBI"] = met.annotation.pop("chebi")
            if "sbo" in met.annotation.keys():
                met.annotation["SBO"] = met.annotation.pop("sbo")
            for annot, val in met.annotation.items():
                if isinstance(val, str):
                    met.annotation[annot] = [val]
    for rxn in model.reactions:
        rxn.name = rxn.name.strip()
        if "sbo" in rxn.annotation.keys():
            rxn.annotation["SBO"] = rxn.annotation.pop("sbo")
        if len(rxn.annotation):
            for annot, val in rxn.annotation.items():
                if isinstance(val, str):
                    rxn.annotation[annot] = [val]
    for gene in model.genes:
        if len(gene.annotation):
            if "ncbigi" in gene.annotation.keys():
                gene.annotation["ncbiprotein"] = gene.annotation.pop("ncbigi")
            for annot, val in gene.annotation.items():
                if isinstance(val, str):
                    gene.annotation[annot] = [val]


@pytest.mark.skipif(scipy is None, reason="scipy unavailable")
@pytest.mark.parametrize(
    "xml_file", ["e_coli_core.xml", "salmonella.xml", "mini_cobra.xml", "mini_fbc2.xml"]
)
# When using a better comparison function, can run test on
# "annotation.xml", "example_notes.xml", "fbc_ex1.xml", "fbc_ex2.xml", "validation.xml"
# "example_notes.xml" contains a group and groups are not yet correctly exported to
# matlab
# "valid_annotation_output.xml" has reaction annotations in a metabolite, so they would
# be thrown out by matlab
def test_compare_xml_to_written_matlab_model(
    data_directory: str, tmpdir: "py.path.local", xml_file: str
) -> None:
    """Verify that xml rewritten as mat file is written and read correctly."""
    xml_model = io.read_sbml_model(join(data_directory, xml_file))
    _fix_xml_annotation_to_identifiers(xml_model)
    mat_output_file = tmpdir.join(xml_file.replace(".xml", ".mat"))
    io.save_matlab_model(
        xml_model, str(mat_output_file)
    )  # lac__D_e_boundary confuses the reading of matlab
    mat_model = io.load_matlab_model(str(mat_output_file))
    assert compare_models(xml_model, mat_model) is None


@pytest.mark.skipif(scipy is None, reason="scipy unavailable")
def test_fail_on_problematic_compartments(data_directory: str) -> None:
    """Test that mat import will fail if there are problems in compartments."""
    with pytest.raises(IOError):
        # AntCore does not have defined compartments
        ant_core_model = io.load_matlab_model(join(data_directory, "AntCore.mat"))
    with pytest.raises(IOError):
        # Ec_iAF1260_flux1 has underscore in compartment names which is not allowed
        Ec_iAF1260_flux1_model = io.load_matlab_model(
            join(data_directory, "Ec_iAF1260_flux1.mat")
        )


@pytest.mark.skipif(scipy is None, reason="scipy unavailable")
def test_mat_model_with_long_compartment_ids(
    data_directory: str, tmpdir: "py.path.local"
) -> None:
    """Test that long compartment IDs like "luSI" are correctly loaded"""
    model_compartments = io.load_matlab_model(join(data_directory, "compartments.mat"))
    assert model_compartments.compartments == {
        "csf": "csf",
        "bcK": "bcK",
        "a": "a",
        "luSI": "luSI",
        "luLI": "luLI",
        "luP": "luP",
        "aL": "aL",
        "fe": "fe",
    }
    assert len(model_compartments.metabolites) == 8
    assert len(model_compartments.reactions) == 15
    for met in model_compartments.metabolites:
        assert met.annotation == {
            "bigg.metabolite": ["glc__D"],
            "cas": ["50-99-7"],
            "kegg.compound": ["C00031"],
            "pubchem.substance": ["3333"],
        }
    output_file = tmpdir.join("compartments.mat")
    io.save_matlab_model(model_compartments, str(output_file))
    model_compartments_reloaded = io.load_matlab_model(str(output_file))
    assert compare_models(model_compartments, model_compartments_reloaded) is None


@pytest.mark.skipif(scipy is None, reason="scipy unavailable")
def test_mat_model_with_no_genes(data_directory: str, tmpdir: "py.path.local") -> None:
    """Test that a model with no genes is loaded and reloaded correctly."""
    model_no_genes = io.load_matlab_model(
        join(data_directory, "cardiac_mit_glcuptake_atpmax.mat")
    )
    assert not len(model_no_genes.genes)
    output_file = tmpdir.join("cardiac_mit_glcuptake_atpmax.mat")
    io.save_matlab_model(model_no_genes, str(output_file))
    model_no_genes_reloaded = io.load_matlab_model(str(output_file))
    assert compare_models(model_no_genes, model_no_genes_reloaded) is None


@pytest.mark.skipif(scipy is None, reason="scipy unavailable")
def test_mat_model_wrong_caps(data_directory: str) -> None:
    """Check that wrong capitalization in matlab field names is processed correctly.

    See https://gist.github.com/akaviaLab/3dcb0eed6563a9d3d1e07198337300ac to create it
    again when needed.
    """
    mat_model = io.load_matlab_model(join(data_directory, "mini.mat"))
    mat_wrong_caps_model = io.load_matlab_model(
        join(data_directory, "mini_wrong_key_caps.mat")
    )
    assert compare_models(mat_model, mat_wrong_caps_model) is None
    assert mat_wrong_caps_model.reactions.get_by_id("LDH_D").annotation == {
        "rhea": ["16369", "16370", "16371", "16372"],
        "metanetx.reaction": ["MNXR101037"],
        "kegg.reaction": ["R00704"],
        "bigg.reaction": ["LDH_D"],
        "ec-code": ["1.1.1.28"],
        "biocyc": ["META:DLACTDEHYDROGNAD-RXN"],
    }
    for rxn in mat_model.reactions.list_attr("id"):
        assert (
            mat_wrong_caps_model.reactions.get_by_id(rxn).annotation
            == mat_model.reactions.get_by_id(rxn).annotation
        )
    assert mat_wrong_caps_model.metabolites.get_by_id("pyr_c").annotation == {
        "seed.compound": ["cpd00020"],
        "unipathway.compound": ["UPC00022"],
        "lipidmaps": ["LMFA01060077"],
        "reactome": ["REACT_113557", "REACT_389680", "REACT_29398"],
        "biocyc": ["PYRUVATE"],
        "pubchem.substance": ["3324"],
        "bigg.metabolite": ["pyr"],
        "cas": ["127-17-3"],
        "hmdb": ["HMDB00243"],
        "kegg.compound": ["C00022"],
    }
    for met in mat_model.metabolites.list_attr("id"):
        assert (
            mat_wrong_caps_model.metabolites.get_by_id(met).annotation
            == mat_model.metabolites.get_by_id(met).annotation
        )
