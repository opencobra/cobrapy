# -*- coding: utf-8 -*-

"""Test functionalities provided by mat.py"""

from __future__ import absolute_import

from builtins import open  # Python 2 unicode compatibility.
from os.path import join
from pickle import load

import pytest
import cobra.io as cio

import helpers


scipy = pytest.importorskip("scipy")


@pytest.fixture(scope="module")
def raven_model(data_directory):
    """Fixture for RAVEN model."""
    with open(join(data_directory, "raven.pickle"),
              "rb", encoding=None) as infile:
        return load(infile)


# @pytest.mark.parametrize("ref_model, filename",
#                          [(pytest.fixture_request("mini_model"),
#                            "mini.mat"),
#                           (pytest.fixture_request("raven_model"),
#                            "raven.mat")])
# TODO: wait for pytest.fixture_request() to get approved
def test_load_matlab_model(data_directory, mini_model, raven_model):
    """Test the reading of MAT model."""
    mini_mat_model = cio.load_matlab_model(join(data_directory, "mini.mat"))
    raven_mat_model = cio.load_matlab_model(join(data_directory, "raven.mat"))
    helpers.assert_equal_models(mini_model, mini_mat_model)
    helpers.assert_equal_models(raven_model, raven_mat_model)


# @pytest.mark.parametrize("model, filename",
#                          [(pytest.fixture_request("mini_model"),
#                            "mini.mat"),
#                           (pytest.fixture_request("raven_model"),
#                            "raven.mat")])
# TODO: wait for pytest.fixture_request() to get approved
def test_save_matlab_model(tmpdir, mini_model, raven_model):
    """Test the writing of MAT model."""
    mini_output_file = tmpdir.join("mini.mat")
    raven_output_file = tmpdir.join("raven.mat")
    # scipy.io.savemat() doesn't support anything other than
    # str or file-stream object, hence the str conversion
    cio.save_matlab_model(mini_model, str(mini_output_file))
    cio.save_matlab_model(raven_model, str(raven_output_file))
    assert mini_output_file.check()
    assert raven_output_file.check()
