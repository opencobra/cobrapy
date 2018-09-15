# -*- coding: utf-8 -*-

"""Test data storage and recovery using pickle."""

from __future__ import absolute_import

from builtins import open  # Python 2 unicode compatibility.
from os.path import join

import pytest

import helpers


PICKLE_MODULES = ["pickle", "cPickle"]


@pytest.mark.parametrize("module", PICKLE_MODULES)
def test_read_pickle(data_directory, mini_model, module):
    """Test the reading of model from pickle."""
    pickle = pytest.importorskip(module)
    with open(join(data_directory, "mini.pickle"),
              "rb", encoding=None) as infile:
        pickle_model = pickle.load(infile)
    helpers.assert_equal_models(mini_model, pickle_model)


@pytest.mark.parametrize("module", PICKLE_MODULES)
def test_write_pickle(tmpdir, mini_model, module):
    """Test the writing of model to pickle."""
    pickle = pytest.importorskip(module)
    output_file = tmpdir.join("mini.pickle")
    with open(str(output_file), "wb", encoding=None) as outfile:
        pickle.dump(mini_model, outfile)
    assert output_file.check()
