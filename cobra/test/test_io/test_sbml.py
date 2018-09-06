# -*- coding: utf-8 -*-

"""Test functionalities provided by sbml.py"""

from __future__ import absolute_import

from os.path import join

import pytest

import cobra.io as cio

import helpers


libsbml = pytest.importorskip("libsbml")


@pytest.fixture(scope="function")
def mini_fbc1_model(data_directory):
    return cio.read_legacy_sbml(join(data_directory, "mini_fbc1.xml"))


@pytest.fixture(scope="function")
def mini_cobra_model(data_directory):
    return cio.read_legacy_sbml(join(data_directory, "mini_cobra.xml"))


# TODO: parametrize the arguments after pytest.fixture_request()
# is approved
def test_read_sbml_model(data_directory, mini_fbc1_model, mini_cobra_model):
    """Test the reading of a model from SBML v2."""
    mini_fbc1 = cio.read_legacy_sbml(join(data_directory, "mini_fbc1.xml"))
    mini_cobra = cio.read_legacy_sbml(join(data_directory, "mini_cobra.xml"))
    helpers.assert_equal_models(mini_fbc1_model, mini_fbc1)
    helpers.assert_equal_models(mini_cobra_model, mini_cobra)


# TODO: parametrize the arguments after pytest.fixture_request()
# is approved
def test_write_sbml_model(tmpdir, mini_fbc1_model, mini_cobra_model):
    """Test the writing of a model to SBML v2."""
    mini_fbc1_output_file = tmpdir.join("mini_fbc1.xml")
    mini_cobra_output_file = tmpdir.join("mini_cobra.xml")

    # convert to str object before passing the filename
    cio.write_legacy_sbml(mini_fbc1_model, str(mini_fbc1_output_file),
                          use_fbc_package=True)
    cio.write_legacy_sbml(mini_cobra_model, str(mini_cobra_output_file),
                          use_fbc_package=False)

    assert mini_fbc1_output_file.check()
    assert mini_cobra_output_file.check()
