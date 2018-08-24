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


@pytest.mark.xfail
@pytest.mark.skipif(libsbml is None, reason="libsbml unavailable.")
@pytest.mark.parametrize("sbml_file", ["mini_fbc1.xml", "mini_cobra.xml"])
def test_read_sbml_model(data_directory, mini_model, sbml_file):
    """Test the reading of a model from SBML2."""
    sbml_model = io.read_legacy_sbml(join(data_directory, sbml_file))
    assert compare_models(mini_model, sbml_model)


@pytest.mark.xfail
@pytest.mark.skipif(libsbml is None, reason="libsbml unavailable.")
@pytest.mark.parametrize("sbml_file, use_fbc", [("mini_fbc1.xml", True),
                                                ("mini_cobra.xml", False)])
def test_write_sbml_model(data_directory, tmpdir, sbml_file, use_fbc):
    """Test the writing of a model to SBML2."""
    input_file = join(data_directory, sbml_file)
    output_file = join(tmpdir, sbml_file)

    sbml_model = io.read_legacy_sbml(input_file)
    io.write_legacy_sbml(sbml_model, output_file, use_fbc_package=use_fbc)
    #TODO: xfail is a temporary workaround until xmldiff is set up.
    assert getsize(input_file) == getsize(output_file)
