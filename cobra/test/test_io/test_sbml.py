# -*- coding: utf-8 -*-

"""Test functionalities provided by sbml.py"""

from __future__ import absolute_import

from os.path import getsize, join

import pytest

from cobra import io
from cobra.test.test_io.conftest import compare_models


try:
    import libsbml
except ImportError:
    libsbml = None


@pytest.fixture(scope="function")
def mini_fbc1_model(data_directory):
    return io.read_legacy_sbml(join(data_directory, "mini_fbc1.xml"))


@pytest.fixture(scope="function")
def mini_cobra_model(data_directory):
    return io.read_legacy_sbml(join(data_directory, "mini_cobra.xml"))


# TODO: parametrize the arguments after pytest.fixture_request()
# is approved
@pytest.mark.skipif(libsbml is None, reason="libsbml unavailable.")
def test_read_sbml_model(data_directory, mini_fbc1_model, mini_cobra_model):
    """Test the reading of a model from SBML v2."""
    mini_fbc1 = io.read_legacy_sbml(join(data_directory, "mini_fbc1.xml"))
    mini_cobra = io.read_legacy_sbml(join(data_directory, "mini_cobra.xml"))
    assert compare_models(mini_fbc1_model, mini_fbc1) is None
    assert compare_models(mini_cobra_model, mini_cobra) is None


# TODO: parametrize the arguments after pytest.fixture_request()
# is approved
@pytest.mark.skipif(libsbml is None, reason="libsbml unavailable.")
def test_write_sbml_model(tmpdir, mini_fbc1_model, mini_cobra_model):
    """Test the writing of a model to SBML v2."""
    mini_fbc1_output_file = tmpdir.join("mini_fbc1.xml")
    mini_cobra_output_file = tmpdir.join("mini_cobra.xml")

    # convert to str object before passing the filename
    io.write_legacy_sbml(mini_fbc1_model, str(mini_fbc1_output_file),
                         use_fbc_package=True)
    io.write_legacy_sbml(mini_cobra_model, str(mini_cobra_output_file),
                         use_fbc_package=False)

    assert mini_fbc1_output_file.check()
    assert mini_cobra_output_file.check()
