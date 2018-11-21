# -*- coding: utf-8 -*-

"""Test functionalities provided by sbml3.py"""

from __future__ import absolute_import

from os.path import join

import pytest
from six import itervalues

from cobra import io
from cobra.test.test_io.conftest import compare_models


@pytest.fixture(scope="function")
def mini_fbc2_model(data_directory):
    """Return mini_fbc2 model."""
    return io.sbml3.read_sbml_model(join(data_directory, "mini_fbc2.xml"))


# Benchmarks
def test_benchmark_read(data_directory, benchmark):
    """Benchmark SBML read."""
    benchmark(io.sbml3.read_sbml_model, join(data_directory, "mini_fbc2.xml"))


def test_benchmark_write(model, benchmark, tmpdir):
    """Benchmark SBML write."""
    benchmark(io.sbml3.write_sbml_model, model, tmpdir.join("-bench"))


# Tests
def test_sbml3_error(data_directory):
    """Test invalid SBML read."""
    filename = join(data_directory, "invalid0.xml")
    with pytest.raises(io.sbml3.CobraSBMLError):
        io.read_sbml_model(filename)


def test_validate_sbml_model(data_directory):
    """Test validation of SBML."""
    # invalid SBML
    for i in range(3):
        filename = join(data_directory, "invalid{}.xml".format(i))
        _, errors = io.sbml3.validate_sbml_model(filename)
        assert all(len(v) >= 1 for v in itervalues(errors)) is False

    # valid SBML
    filename = join(data_directory, "mini_fbc2.xml")
    _, errors = io.sbml3.validate_sbml_model(filename)
    assert all(len(v) == 0 for v in itervalues(errors))


@pytest.mark.parametrize("sbml_file", ["mini_fbc2.xml", "mini_fbc2.xml.gz",
                                       "mini_fbc2.xml.bz2"])
def test_read_sbml_model(data_directory, mini_model, sbml_file):
    """Test the reading of a model from SBML3."""
    sbml3_model = io.read_sbml_model(join(data_directory, sbml_file))
    assert compare_models(mini_model, sbml3_model) is None


@pytest.mark.parametrize("ext", [".xml", ".xml.gz", ".xml.bz2"])
def test_write_sbml_model(tmpdir, mini_fbc2_model, ext):
    """Test the writing of model to SBML3."""
    output_file = tmpdir.join("mini_fbc2{}".format(ext))
    io.write_sbml_model(mini_fbc2_model, output_file)
    assert output_file.check()
