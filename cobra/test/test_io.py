# -*- coding: utf-8 -*-
from __future__ import absolute_import

from collections import namedtuple
from functools import partial
from os import unlink
from os.path import join, split
from pickle import dump, load
from tempfile import gettempdir
from warnings import warn

import pytest
from six import iteritems

from cobra import io


def write_legacy_sbml_placeholder():
    pass


try:
    import scipy
except ImportError:
    scipy = None
try:
    import libsbml

    write_legacy_sbml = io.write_legacy_sbml
except ImportError:
    libsbml = None
    write_legacy_sbml = write_legacy_sbml_placeholder
try:
    import jsonschema
except ImportError:
    jsonschema = None
try:
    import cPickle

    cload = cPickle.load
    cdump = cPickle.dump
except ImportError:
    cPickle = None
    cload = None
    cdump = None


def validate_json(filename):
    with open(filename, "r") as infile:
        loaded = io.json.json.load(infile)
    if jsonschema is None:
        warn("jsonschema not installed")
    else:
        jsonschema.validate(loaded, io.json.json_schema)


def read_pickle(filename, load_function=load):
    with open(filename, "rb") as infile:
        return load_function(infile)


def write_pickle(model, filename, dump_function=dump):
    with open(filename, "wb") as outfile:
        dump_function(model, outfile)


IOTrial = namedtuple('IOTrial',
                     ['name', 'reference_file', 'test_file', 'read_function',
                      'write_function', 'validation_function'])
trials = [IOTrial('fbc2', 'mini.pickle', 'mini_fbc2.xml',
                  io.read_sbml_model, io.write_sbml_model,
                  io.sbml3.validate_sbml_model),
          IOTrial('fbc2Gz', 'mini.pickle', 'mini_fbc2.xml.gz',
                  io.read_sbml_model, io.write_sbml_model, None),
          IOTrial('fbc2Bz2', 'mini.pickle', 'mini_fbc2.xml.bz2',
                  io.read_sbml_model, io.write_sbml_model, None),
          pytest.mark.skipif("not libsbml")(
              IOTrial('fbc1', 'mini.pickle', 'mini_fbc1.xml',
                      io.read_sbml_model,
                      partial(write_legacy_sbml, use_fbc_package=True), None)),
          pytest.mark.skipif("not libsbml")(
              IOTrial('cobra', 'mini.pickle', 'mini_cobra.xml',
                      io.read_sbml_model,
                      partial(write_legacy_sbml, use_fbc_package=False),
                      None)),
          pytest.mark.skipif("not scipy")(
              IOTrial('mat', 'mini.pickle', 'mini.mat',
                      io.load_matlab_model, io.save_matlab_model, None)),
          pytest.mark.skipif("not scipy")(
              IOTrial('raven-mat', 'raven.pickle', 'raven.mat',
                      io.load_matlab_model, io.save_matlab_model, None)),
          IOTrial('json', 'mini.pickle', 'mini.json',
                  io.load_json_model, io.save_json_model, validate_json),
          IOTrial('yaml', 'mini.pickle', 'mini.yml',
                  io.load_yaml_model, io.save_yaml_model, None),
          IOTrial('pickle', 'mini.pickle', 'mini.pickle',
                  read_pickle, write_pickle, None),
          pytest.mark.skipif("not cPickle")(
              IOTrial('pickle', 'mini.pickle', 'mini.pickle',
                      partial(read_pickle, load_function=cload),
                      partial(write_pickle, dump_function=cdump), None))
          ]
trial_names = [node.name for node in trials]


@pytest.mark.skipif(scipy is not None, reason='scipy available')
def raise_scipy_errors():
    with pytest.raises(ImportError):
        io.save_matlab_model(None, 'test')
    with pytest.raises(ImportError):
        io.load_matlab_model('test')


@pytest.mark.skipif(libsbml is not None, reason='libsbml available')
def raise_libsbml_errors():
    with pytest.raises(ImportError):
        io.read_sbml_model('test')
    with pytest.raises(ImportError):
        io.write_sbml_model(None, 'test')
    with pytest.raises(ImportError):
        io.load_matlab_model('test')
    with pytest.raises(ImportError):
        io.write_legacy_sbml(None, 'test')
    with pytest.raises(ImportError):
        io.read_legacy_sbml(None, 'test')


@pytest.fixture(scope="module", params=trials, ids=trial_names)
def io_trial(request, data_directory):
    with open(join(data_directory, request.param.reference_file),
              "rb") as infile:
        reference_model = load(infile)
    test_model = request.param.read_function(join(data_directory,
                                                  request.param.test_file))
    test_output_filename = join(gettempdir(),
                                split(request.param.test_file)[-1])
    # test writing the model within a context with a non-empty stack
    with test_model:
        test_model.objective = test_model.objective
        request.param.write_function(test_model, test_output_filename)
    reread_model = request.param.read_function(test_output_filename)
    unlink(test_output_filename)
    return request.param.name, reference_model, test_model, reread_model


class TestCobraIO:
    @classmethod
    def compare_models(cls, name, model1, model2):
        assert len(model1.reactions) == len(model2.reactions)
        assert len(model1.metabolites) == len(model2.metabolites)
        assert model1.objective.direction == model2.objective.direction
        for attr in ("id", "name", "lower_bound", "upper_bound",
                     "objective_coefficient", "gene_reaction_rule"):
            assert getattr(model1.reactions[0], attr) == getattr(
                model2.reactions[0], attr)
            assert getattr(model1.reactions[5], attr) == getattr(
                model2.reactions[5], attr)
            assert getattr(model1.reactions[-1], attr) == getattr(
                model2.reactions[-1], attr)
        for attr in ("id", "name", "compartment", "formula", "charge"):
            assert getattr(model1.metabolites[0], attr) == getattr(
                model2.metabolites[0], attr)
            assert getattr(model1.metabolites[5], attr) == getattr(
                model2.metabolites[5], attr)
            assert getattr(model1.metabolites[-1], attr) == getattr(
                model2.metabolites[-1], attr)
        assert len(model1.reactions[0].metabolites) == len(
            model2.reactions[0].metabolites)
        assert len(model1.reactions[8].metabolites) == len(
            model2.reactions[8].metabolites)
        assert len(model1.reactions[-1].metabolites) == len(
            model2.reactions[-1].metabolites)
        assert len(model1.genes) == len(model2.genes)
        # ensure they have the same solution max
        solution1 = model1.optimize()
        solution2 = model2.optimize()
        assert abs(solution1.f - solution2.f) < 0.001
        # ensure the references are correct
        assert model2.metabolites[0]._model is model2
        assert model2.reactions[0]._model is model2
        assert model2.genes[0]._model is model2

    @classmethod
    def extra_comparisons(cls, name, model1, model2):
        assert model1.compartments == model2.compartments
        assert dict(model1.metabolites[4].annotation) == dict(
            model2.metabolites[4].annotation)
        assert dict(model1.reactions[4].annotation) == dict(
            model2.reactions[4].annotation)
        assert dict(model1.genes[5].annotation) == dict(
            model2.genes[5].annotation)
        for attr in ("id", "name"):
            assert getattr(model1.genes[0], attr) == getattr(model2.genes[0],
                                                             attr)
            assert getattr(model1.genes[10], attr) == getattr(model2.genes[10],
                                                              attr)
            assert getattr(model1.genes[-1], attr) == getattr(model2.genes[-1],
                                                              attr)

    def test_read_1(self, io_trial):
        name, reference_model, test_model, _ = io_trial
        if name in ['fbc1']:
            pytest.xfail('not supported')
        self.compare_models(name, reference_model, test_model)

    def test_read_2(self, io_trial):
        name, reference_model, test_model, _ = io_trial
        if name in ['fbc1', 'mat', 'cobra', 'raven-mat']:
            pytest.xfail('not supported')
        self.extra_comparisons(name, reference_model, test_model)

    def test_write_1(self, io_trial):
        name, _, test_model, reread_model = io_trial
        if name in ['fbc1', 'raven-mat']:
            pytest.xfail('not supported')
        self.compare_models(name, test_model, reread_model)

    def test_write_2(self, io_trial):
        name, _, test_model, reread_model = io_trial
        if name in ['fbc1', 'mat', 'cobra', 'raven-mat']:
            pytest.xfail('not supported')
        self.extra_comparisons(name, test_model, reread_model)


def test_benchmark_read(data_directory, benchmark):
    benchmark(io.sbml3.read_sbml_model, join(data_directory, 'mini_fbc2.xml'))


def test_benchmark_write(model, benchmark):
    benchmark(io.sbml3.write_sbml_model, model, join(gettempdir(), "-bench"))


@pytest.mark.parametrize("trial", trials)
def test_validate(trial, data_directory):
    if trial.validation_function is None:
        pytest.skip('not implemented')
    test_file = join(data_directory, trial.test_file)
    trial.validation_function(test_file)


@pytest.mark.parametrize("trial", trials)
def test_read_nonexistent(trial):
    pytest.raises(IOError, trial.read_function, "fake_file")


def test_sbml_error(data_directory):
    filename = join(data_directory, "invalid0.xml")
    with pytest.raises(io.sbml3.CobraSBMLError):
        io.read_sbml_model(filename)


def test_bad_validation(data_directory):
    for i in range(3):
        filename = join(data_directory, "invalid%d.xml" % i)
        m, errors = io.sbml3.validate_sbml_model(filename)
        assert any(len(v) >= 1 for v in iteritems(errors))
