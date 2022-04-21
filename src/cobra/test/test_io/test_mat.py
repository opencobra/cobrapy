"""Test functionalities of I/O in MATLAB (.mat) format."""

import pathlib
from os.path import join
from pickle import load
from typing import TYPE_CHECKING

import pytest

from cobra import io
from cobra.test.test_io.conftest import compare_models


try:
    import scipy
except ImportError:
    scipy = None


if TYPE_CHECKING:
    import py.path

    from cobra import Model


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
    data_directory: str, mini_model: "Model", raven_model: "Model"
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
def test_read_rewrite_matlab_model(tmpdir: "py.path.local", data_directory: str) -> None:
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


@pytest.mark.skipif(scipy is None, reason="scipy unavailable")
def test_compare_xml_to_written_matlab_model(data_directory: str, tmpdir: "py.path.local") -> None:
    """Verify that xml rewritten as mat file is written and read correctly."""
    for xml_file in pathlib.Path(data_directory).glob('.xml'):
        xml_model = io.read_sbml_model(str(xml_file))
        mat_output_file = tmpdir.join(xml_file.name.replace('.xml', '.mat'))
        io.save_matlab_model(xml_model, str(mat_output_file))
        mat_model = io.load_matlab_model(str(mat_output_file))
        assert compare_models(xml_model, mat_model) is None


@pytest.mark.skipif(scipy is None, reason="scipy unavailable")
def test_fail_on_problematic_compartments(data_directory: str) -> None:
    """Test that mat import will fail if there are problems in compartments."""
    with pytest.raises(IOError):
        # AntCore does not have defined compartments
        ant_core_model = io.load_matlab_model(join(data_directory, 'AntCore.mat'))
    with pytest.raises(IOError):
        # Ec_iAF1260_flux1 has underscore in compartment names which is not allowed
        Ec_iAF1260_flux1_model = io.load_matlab_model(join(data_directory,
                                                           'Ec_iAF1260_flux1.mat'))


# Test for lots of letters in compartments (short Harvey)

# Test for lots of annotations. short RECON3?

@pytest.mark.skipif(scipy is None, reason="scipy unavailable")
def test_mat_model_with_no_genes(data_directory: str, tmpdir: "py.path.local") -> None:
    """Test that a model with no genes is loaded and reloaded correctly."""
    model_no_genes = io.load_matlab_model(join(data_directory,
                                               'cardiac_mit_glcuptake_atpmax.mat'))
    assert not len(model_no_genes.genes)
    output_file = tmpdir.join("cardiac_mit_glcuptake_atpmax.mat")
    io.save_matlab_model(model_no_genes, str(output_file))
    model_no_genes_reloaded = io.load_matlab_model(str(output_file))
    assert compare_models(model_no_genes, model_no_genes_reloaded) is None
