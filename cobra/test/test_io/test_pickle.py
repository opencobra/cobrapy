# -*- coding: utf-8 -*-

"""Test data storage and recovery using pickle."""

from __future__ import absolute_import

from os.path import join
from pickle import dump, load

import pytest

from cobra.test.test_io.conftest import compare_models


try:
    import cPickle
    cload = cPickle.load
    cdump = cPickle.dump
except ImportError:
    cload = None
    cdump = None


@pytest.mark.parametrize("load_function", [load, cload])
def test_read_pickle(data_directory, mini_model, load_function):
    """Test the reading of model from pickle."""
    if load_function is None:
        pytest.skip()

    with open(join(data_directory, "mini.pickle"), "rb") as infile:
        pickle_model = load_function(infile)

    assert compare_models(mini_model, pickle_model) is None


@pytest.mark.parametrize("dump_function", [dump, cdump])
def test_write_pickle(tmpdir, mini_model, dump_function):
    """Test the writing of model to pickle."""
    if dump_function is None:
        pytest.skip()

    output_file = tmpdir.join("mini.pickle")
    with open(str(output_file), "wb") as outfile:
        dump_function(mini_model, outfile)

    assert output_file.check()
