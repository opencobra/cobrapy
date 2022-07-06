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
@pytest.mark.parametrize(
    "fix_ref_model, fix_directory, filename",
    [
        ("mini_model", "cobra_data_directory", "mini.mat"),
        ("raven_model", "data_directory", "raven.mat"),
    ],
)
def test_load_matlab_model(
    compare_models: Callable,
    fix_ref_model: str,
    fix_directory: str,
    filename: str,
    request: pytest.FixtureRequest,
) -> None:
    """Test the reading of MAT model.

    Will check Path and str for each model.

    Parameters
    ----------
    compare_models : Callable
        A callable function to compare models.
    fix_ref_model: str
        Name of reference model fixture which will be requested.
    fix_directory: str
        Name of directory fixture which will be requested.
    filename: str
        Filename to use as a parameter
    request: FixtureRequest
        Will be used to request fixtures.

    """
    current_model = load_matlab_model(
        (request.getfixturevalue(fix_directory) / filename).resolve()
    )
    assert compare_models(request.getfixturevalue(fix_ref_model), current_model) is None
    current_model = load_matlab_model(
        str((request.getfixturevalue(fix_directory) / filename).resolve())
    )
    assert compare_models(request.getfixturevalue(fix_ref_model), current_model) is None
    with request.getfixturevalue(fix_directory).joinpath(filename).open("rb") as file_h:
        current_model = load_matlab_model(file_h)
        assert (
            compare_models(request.getfixturevalue(fix_ref_model), current_model)
            is None
        )


# @pytest.mark.xfail(reason="localPath not supported yet")
@pytest.mark.skipif(scipy is None, reason="scipy unavailable")
@pytest.mark.parametrize(
    "fix_ref_model, filename",
    [("mini_model", "mini.mat"), ("raven_model", "raven.mat")],
)
def test_save_matlab_model(tmp_path: Path, fix_ref_model, filename, request) -> None:
    """Test the writing of MAT model.

    Parameters
    ----------
    tmp_path : pathlib.Path
        The path to the temporary test assets store.
    fix_ref_model: str
        Name of reference model fixture which will be requested.
    filename: str
        Filename to use as a parameter
    request: FixtureRequest
        Will be used to request fixtures.
    """
    output_file = tmp_path / filename
    save_matlab_model(
        request.getfixturevalue(fix_ref_model), str(output_file.resolve())
    )
    assert output_file.exists()
    output_file.unlink()
    save_matlab_model(request.getfixturevalue(fix_ref_model), output_file.resolve())
    assert output_file.exists()
    output_file.unlink()
    with output_file.open("wb") as file_h:
        save_matlab_model(request.getfixturevalue(fix_ref_model), file_h)
    assert output_file.exists()


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
    save_matlab_model(model, filepath.resolve())
    read = load_matlab_model(filepath.resolve())
    assert read.reactions[0].bounds == (-1e6, 1e6)


@pytest.mark.skipif(scipy is None, reason="scipy unavailable")
@pytest.mark.parametrize(
    "fix_directory, filename",
    [("cobra_data_directory", "mini.mat"), ("data_directory", "raven.mat")],
)
def test_read_rewrite_matlab_model(
    compare_models: Callable,
    tmp_path: Path,
    fix_directory: str,
    filename: str,
    request: pytest.FixtureRequest,
) -> None:
    """Verify that rewritten matlab model is identical to original.

    Parameters
    ----------
    compare_models : Callable
        A callable to compare models.
    tmp_path : pathlib.Path
        The path to the temporary test assets store.
    fix_directory: str
        Name of directory fixture which will be requested.
    filename: str
        Filename to use as a parameter
    request: FixtureRequest
        Will be used to request fixtures.

    """
    current_model = load_matlab_model(request.getfixturevalue(fix_directory) / filename)
    output_file = tmp_path.joinpath(filename)
    save_matlab_model(current_model, output_file)
    mat_model_reload = load_matlab_model(output_file)
    assert compare_models(current_model, mat_model_reload) is None


@pytest.mark.skipif(scipy is None, reason="scipy unavailable")
@pytest.mark.parametrize(
    "dirname, xml_file",
    [
        ("data_directory", "e_coli_core.xml"),
        ("cobra_data_directory", "salmonella.xml.gz"),
        ("cobra_data_directory", "mini_cobra.xml"),
        ("data_directory", "mini_fbc2.xml"),
    ],
)
# When using a better comparison function, can run test on
# "annotation.xml", "example_notes.xml", "fbc_ex1.xml", "fbc_ex2.xml", "validation.xml"
# "valid_annotation_output.xml" has reaction annotations in a metabolite, so they would
# be thrown out by matlab
def test_compare_xml_to_written_matlab_model(
    compare_models: Callable,
    tmp_path: Path,
    dirname: str,
    xml_file: str,
    request: pytest.FixtureRequest,
) -> None:
    """Verify that xml rewritten as mat file is written and read correctly.

    Parameters
    ----------
    compare_models : Callable
        A callable to compare models.
    tmp_path : pathlib.Path
        The path to the temporary test assets store.
    xml_file : str
        The name of the XML file to compare against.

    """
    xml_model = read_sbml_model(request.getfixturevalue(dirname) / xml_file)
    mat_output_file = tmp_path / xml_file.replace(".xml", ".mat")
    save_matlab_model(xml_model, mat_output_file)
    mat_model = load_matlab_model(mat_output_file)
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
        load_matlab_model(data_directory / "AntCore.mat")
    with pytest.raises(IOError):
        # Ec_iAF1260_flux1 has underscore in compartment names which is not allowed
        load_matlab_model(data_directory / "Ec_iAF1260_flux1.mat")


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
    model_compartments = load_matlab_model(data_directory / "compartments.mat")
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
    save_matlab_model(model_compartments, output_file.resolve())
    model_compartments_reloaded = load_matlab_model(output_file.resolve())
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
        data_directory / "cardiac_mit_glcuptake_atpmax.mat"
    )
    assert not len(model_no_genes.genes)
    output_file = tmp_path / "cardiac_mit_glcuptake_atpmax.mat"
    save_matlab_model(model_no_genes, output_file.resolve())
    model_no_genes_reloaded = load_matlab_model(output_file.resolve())
    assert compare_models(model_no_genes, model_no_genes_reloaded) is None


@pytest.mark.skipif(scipy is None, reason="scipy unavailable")
def test_mat_model_wrong_caps(
    compare_models: Callable, data_directory: Path, cobra_data_directory: Path
) -> None:
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
    mat_model = load_matlab_model((cobra_data_directory / "mini.mat").resolve())
    mat_wrong_caps_model = load_matlab_model(
        (data_directory / "mini_wrong_key_caps.mat").resolve()
    )
    assert compare_models(mat_model, mat_wrong_caps_model) is None
    assert mat_wrong_caps_model.reactions.get_by_id("LDH_D").annotation == {
        "rhea": ["16369", "16370", "16371", "16372"],
        "metanetx.reaction": ["MNXR101037"],
        "kegg.reaction": ["R00704"],
        "bigg.reaction": ["LDH_D"],
        "ec-code": ["1.1.1.28"],
        "biocyc": ["META:DLACTDEHYDROGNAD-RXN"],
        "sbo": ["SBO:0000375"],
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
        "chebi": [
            "CHEBI:15361",
            "CHEBI:14987",
            "CHEBI:8685",
            "CHEBI:32816",
            "CHEBI:45253",
            "CHEBI:26466",
            "CHEBI:26462",
        ],
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
