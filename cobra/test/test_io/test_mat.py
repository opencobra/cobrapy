# -*- coding: utf-8 -*-

"""Test functionalities provided by mat.py"""

from __future__ import absolute_import

from os.path import getsize, join
from pickle import load

import pytest
from cobra import io
from cobra.test.test_io.conftest import compare_models

try:
    import scipy
except ImportError:
    scipy = None


@pytest.fixture(scope="function")
def raven_model(data_directory):
    """Fixture for RAVEN model."""
    with open(join(data_directory, "raven.pickle"), "rb") as infile:
        return load(infile)


@pytest.mark.xfail
@pytest.mark.skipif(scipy is None, reason='SciPy unavailable')
# @pytest.mark.parametrize("ref_model, filename",
#                          [(pytest.fixture_request("mini_model"),
#                            "mini.mat"),
#                           (pytest.fixture_request("raven_model"),
#                            "raven.mat")])
# waiting for pytest.fixture_request() to get approved
def test_load_matlab_model(data_directory, mini_model, raven_model):
    """Test the reading of MAT model."""
    mini_mat_model = io.load_matlab_model(join(data_directory, "mini.mat"))
    raven_mat_model = io.load_matlab_model(join(data_directory, "raven.mat"))
    #TODO: fix genes ordering from .mat files, till then mark xfail
    assert compare_models(mini_model, mini_mat_model) is None
    assert compare_models(raven_model, raven_mat_model) is None


@pytest.mark.skipif(scipy is None, reason='SciPy unavailable')
@pytest.mark.parametrize("mat_file", ["mini.mat", "raven.mat"])
def test_save_matlab_model(data_directory, tmpdir, mat_file):
    """Test the writing of MAT model."""
    if mat_file == "raven.mat":
        pytest.skip()
    input_file = join(data_directory, mat_file)
    output_file = join(tmpdir, mat_file)

    mat_model = io.load_matlab_model(input_file)
    io.save_matlab_model(mat_model, output_file)
    #TODO: fix getsize for raven, till then skip
    assert getsize(input_file) == getsize(output_file)
