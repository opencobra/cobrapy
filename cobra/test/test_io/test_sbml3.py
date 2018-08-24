# -*- coding: utf-8 -*-

"""Test functionalities provided by sbml3.py"""

from __future__ import absolute_import

from os.path import getsize, join

from six import itervalues

import pytest
from cobra import io
from cobra.test.test_io.conftest import compare_models


# Benchmarks
def test_benchmark_read(data_directory, benchmark):
    """Benchmark SBML read."""
    benchmark(io.sbml3.read_sbml_model, join(data_directory, 'mini_fbc2.xml'))


def test_benchmark_write(model, benchmark, tmpdir):
    """Benchmark SBML write."""
    benchmark(io.sbml3.write_sbml_model, model, join(tmpdir, "-bench"))


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


@pytest.mark.xfail
@pytest.mark.parametrize("sbml_file", ["mini_fbc2.xml", "mini_fbc2.xml.gz",
                                       "mini_fbc2.xml.bz2"])
def test_write_sbml_model(data_directory, tmpdir, sbml_file):
    """Test the writing of model to SBML3."""
    input_file = join(data_directory, sbml_file)
    output_file = join(tmpdir, sbml_file)

    sbml3_model = io.read_sbml_model(input_file)
    io.write_sbml_model(sbml3_model, output_file)
    #TODO: mark xfail until xml diff tool is ready
    assert getsize(input_file) == getsize(output_file)
