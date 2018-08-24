# -*- coding: utf-8 -*-

"""Test data storage and recovery using pickle."""

from __future__ import absolute_import

from os.path import getsize, join
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


@pytest.mark.xfail
@pytest.mark.parametrize("dump_function", [dump, cdump])
def test_write_pickle(data_directory, tmpdir, dump_function):
    """Test the writing of model to pickle."""
    if dump_function is None:
        pytest.skip()

    input_file = join(data_directory, "mini.pickle")
    output_file = join(tmpdir, "mini.pickle")

    # reading
    with open(input_file, "rb") as infile:
        pickle_model = load(infile)

    # writing
    with open(output_file, "wb") as outfile:
        dump(pickle_model, outfile)

    assert getsize(input_file) == getsize(output_file)
