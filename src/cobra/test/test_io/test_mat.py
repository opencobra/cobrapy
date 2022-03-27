"""Test functionalities of I/O in MATLAB (.mat) format."""

from os.path import join
from pickle import load
from typing import TYPE_CHECKING

import pathlib
import pytest

from cobra import io
from cobra.test.test_io.conftest import compare_models


try:
    import scipy
except ImportError:
    scipy = None


if TYPE_CHECKING:
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
    tmpdir: pathlib.Path, mini_model: "Model", raven_model: "Model"
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
def test_large_bounds(tmpdir: pathlib.Path, model: "Model"):
    """Verify that mat bounds don't get broken by the config defaults."""
    model.reactions[0].bounds = -1e6, 1e6
    filepath = str(tmpdir.join("model.mat"))
    io.save_matlab_model(model, filepath)
    read = io.load_matlab_model(filepath)
    assert read.reactions[0].bounds == (-1e6, 1e6)
