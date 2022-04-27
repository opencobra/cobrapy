"""Test functionalities of I/O in MATLAB (.mat) format."""

from pathlib import Path
from pickle import load
from typing import TYPE_CHECKING, Callable

import pytest

from cobra.io import load_matlab_model, read_sbml_model, save_matlab_model


try:
    import scipy
except ImportError:
    scipy = None


if TYPE_CHECKING:
    from cobra import Model


@pytest.fixture(scope="function")
def raven_model(data_directory: Path) -> "Model":
    """Fixture for RAVEN model."""
    with open(data_directory / "raven.pickle", "rb") as infile:
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
    data_directory: Path,
    mini_model: "Model",
    raven_model: "Model",
) -> None:
    """Test the reading of MAT model.

    Parameters
    ----------
    compare_models : Callable
        A callable to compare models.
    data_directory : pathlib.Path
        The path to the test data directory.

    """
    mini_mat_model = load_matlab_model(str((data_directory / "mini.mat").resolve()))
    raven_mat_model = load_matlab_model(str((data_directory / "raven.mat").resolve()))
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
    tmp_path: Path, mini_model: "Model", raven_model: "Model"
) -> None:
    """Test the writing of MAT model.

    Parameters
    ----------
    tmp_path : pathlib.Path
        The path to the temporary test assets store.
    mini_model : cobra.Model
        The mini model.
    raven_model : cobra.Model
        The RAVEN model.

    """
    mini_output_file = tmp_path / "mini.mat"
    raven_output_file = tmp_path / "raven.mat"
    # scipy.io.savemat() doesn't support anything other than
    # str or file-stream object, hence the str conversion
    save_matlab_model(mini_model, str(mini_output_file.resolve()))
    save_matlab_model(raven_model, str(raven_output_file.resolve()))
    assert mini_output_file.exists()
    assert raven_output_file.exists()


@pytest.mark.skipif(scipy is None, reason="scipy unavailable")
def test_large_bounds(tmp_path: Path, model: "Model") -> None:
    """Verify that mat bounds don't get broken by the config defaults.

    Parameters
    ----------
    tmp_path : pathlib.Path
        The path to the temporary test assets store.
    model : cobra.Model
        The "textbook" model.

    """
    model.reactions[0].bounds = -1e6, 1e6
    filepath = tmp_path / "model.mat"
    save_matlab_model(model, str(filepath.resolve()))
    read = load_matlab_model(str(filepath.resolve()))
    assert read.reactions[0].bounds == (-1e6, 1e6)


@pytest.mark.skipif(scipy is None, reason="scipy unavailable")
def test_read_rewrite_matlab_model(
    compare_models: Callable, tmp_path: Path, data_directory: Path
) -> None:
    """Verify that rewritten matlab model is identical to original.

    Parameters
    ----------
    compare_models : Callable
        A callable to compare models.
    tmp_path : pathlib.Path
        The path to the temporary test assets store.
    data_directory : pathlib.Path
        The path to the test data directory.

    """
    mini_mat_model = load_matlab_model(str((data_directory / "mini.mat").resolve()))
    raven_mat_model = load_matlab_model(str((data_directory / "raven.mat").resolve()))
    mini_output_file = tmp_path.joinpath("mini.mat")
    raven_output_file = tmp_path.joinpath("raven.mat")
    # scipy.io.savemat() doesn't support anything other than
    # str or file-stream object, hence the str conversion
    save_matlab_model(mini_mat_model, str(mini_output_file))
    save_matlab_model(raven_mat_model, str(raven_output_file))
    mini_mat_model_reload = load_matlab_model(str(mini_output_file))
    raven_mat_model_reload = load_matlab_model(str(raven_output_file))
    assert compare_models(mini_mat_model, mini_mat_model_reload) is None
    assert compare_models(raven_mat_model, raven_mat_model_reload) is None


def _fix_xml_annotation_to_identifiers(model: "Model") -> None:
    """Fix XML annotations to respect identifiers.org .

    This function will fix the dict keys of annotations to match identifiers.org.
    Eventually, the XML models should be fixed and cobrapy should be strict, but this is
    part of SBML rewriting of annotations
    see: https://github.com/opencobra/cobrapy/issues/684

    It also changes met formulas from empty string to None (which is the default
    when creating a metabolite with no fomula given) and strips spaces from reaction
    names.

    Parameters
    ----------
    model : cobra.Model
        The model to fix.

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
    compare_models: Callable,
    data_directory: Path,
    tmp_path: Path,
    xml_file: str,
) -> None:
    """Verify that xml rewritten as mat file is written and read correctly.

    Parameters
    ----------
    compare_models : Callable
        A callable to compare models.
    data_directory : pathlib.Path
        The path to the test data directory.
    tmp_path : pathlib.Path
        The path to the temporary test assets store.
    xml_file : str
        The name of the XML file to compare against.

    """
    xml_model = read_sbml_model(str((data_directory / xml_file).resolve()))
    mat_output_file = tmp_path / xml_file.replace(".xml", ".mat")
    save_matlab_model(
        xml_model, str(mat_output_file.resolve())
    )  # lac__D_e_boundary confuses the reading of matlab
    mat_model = load_matlab_model(str(mat_output_file.resolve()))
    assert compare_models(xml_model, mat_model) is None


@pytest.mark.skipif(scipy is None, reason="scipy unavailable")
def test_fail_on_problematic_compartments(data_directory: Path) -> None:
    """Test that mat import will fail if there are problems in compartments.

    Parameters
    ----------
    data_directory : pathlib.Path
        The path to the test data directory.

    """
    with pytest.raises(IOError):
        # AntCore does not have defined compartments
        load_matlab_model(str((data_directory / "AntCore.mat").resolve()))
    with pytest.raises(IOError):
        # Ec_iAF1260_flux1 has underscore in compartment names which is not allowed
        load_matlab_model(str((data_directory / "Ec_iAF1260_flux1.mat").resolve()))


@pytest.mark.skipif(scipy is None, reason="scipy unavailable")
def test_mat_model_with_long_compartment_ids(
    compare_models: Callable, data_directory: Path, tmp_path: Path
) -> None:
    """Test that long compartment IDs like "luSI" are correctly loaded.

    Parameters
    ----------
    compare_models : Callable
        A callable to compare models.
    data_directory : pathlib.Path
        The path to the test data directory.
    tmp_path : pathlib.Path
        The path to the temporary test assets store.

    """
    model_compartments = load_matlab_model(
        str((data_directory / "compartments.mat").resolve())
    )
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
    output_file = tmp_path / "compartments.mat"
    save_matlab_model(model_compartments, str(output_file.resolve()))
    model_compartments_reloaded = load_matlab_model(str(output_file.resolve()))
    assert compare_models(model_compartments, model_compartments_reloaded) is None


@pytest.mark.skipif(scipy is None, reason="scipy unavailable")
def test_mat_model_with_no_genes(
    compare_models: Callable, data_directory: Path, tmp_path: Path
) -> None:
    """Test that a model with no genes is loaded and reloaded correctly.

    Parameters
    ----------
    compare_models : Callable
        A callable to compare models.
    data_directory : pathlib.Path
        The path to the test data directory.
    tmp_path : pathlib.Path
        The path to the temporary test assets store.

    """
    model_no_genes = load_matlab_model(
        str((data_directory / "cardiac_mit_glcuptake_atpmax.mat").resolve())
    )
    assert not len(model_no_genes.genes)
    output_file = tmp_path / "cardiac_mit_glcuptake_atpmax.mat"
    save_matlab_model(model_no_genes, str(output_file.resolve()))
    model_no_genes_reloaded = load_matlab_model(str(output_file.resolve()))
    assert compare_models(model_no_genes, model_no_genes_reloaded) is None


@pytest.mark.skipif(scipy is None, reason="scipy unavailable")
def test_mat_model_wrong_caps(compare_models: Callable, data_directory: Path) -> None:
    """Check that wrong capitalization in matlab field names is processed correctly.

    See https://gist.github.com/akaviaLab/3dcb0eed6563a9d3d1e07198337300ac to create it
    again when needed.

    Parameters
    ----------
    compare_models : Callable
        A callable to compare models.
    data_directory : pathlib.Path
        The path to the test data directory.

    """
    mat_model = load_matlab_model(str(Path(data_directory / "mini.mat").resolve()))
    mat_wrong_caps_model = load_matlab_model(
        str(Path(data_directory, "mini_wrong_key_caps.mat").resolve())
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
