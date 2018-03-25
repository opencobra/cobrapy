# -*- coding: utf-8 -*-
"""
Testing SBML functionality based on libsbml.
"""

from __future__ import absolute_import


from os import unlink
from os.path import join, split
from pickle import load
from tempfile import gettempdir
from collections import namedtuple
from functools import partial
from warnings import warn
from six import iteritems


import pytest
from cobra.io import sbmlnew

try:
    import libsbml
except ImportError:
    libsbml = None

try:
    import jsonschema
except ImportError:
    jsonschema = None

# ----------------------------------
# Definition of SBML files to test
# ----------------------------------
IOTrial = namedtuple('IOTrial',
                     ['name', 'reference_file', 'test_file', 'read_function',
                      'write_function', 'validation_function'])
trials = [IOTrial('fbc2', 'mini.pickle', 'mini_fbc2.xml',
                  sbmlnew.read_sbml_model, sbmlnew.write_sbml_model,
                  sbmlnew.validate_sbml_model),
          IOTrial('fbc2Gz', 'mini.pickle', 'mini_fbc2.xml.gz',
                  sbmlnew.read_sbml_model, sbmlnew.write_sbml_model, None),
          IOTrial('fbc2Bz2', 'mini.pickle', 'mini_fbc2.xml.bz2',
                  sbmlnew.read_sbml_model, sbmlnew.write_sbml_model, None),
          IOTrial('fbc1', 'mini.pickle', 'mini_fbc1.xml',
                  sbmlnew.read_sbml_model, sbmlnew.write_sbml_model, None),
          IOTrial('cobra', None, 'mini_cobra.xml',
                  sbmlnew.read_sbml_model, sbmlnew.write_sbml_model, None),
          ]
trial_names = [node.name for node in trials]


@pytest.mark.parametrize("trial", trials)
def test_validate(trial, data_directory):
    """ Test validation function. """
    if trial.validation_function is None:
        pytest.skip('not implemented')
    test_file = join(data_directory, trial.test_file)
    trial.validation_function(test_file)


class TestCobraIO:
    """ Tests the read and write functions. """

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

        print(model2.genes)
        print(model2.genes[0], print(model2.genes[0]._model))
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
        if reference_model:
            self.compare_models(name, reference_model, test_model)

    def test_read_2(self, io_trial):
        name, reference_model, test_model, _ = io_trial
        if name in ['fbc1', 'mat', 'cobra', 'raven-mat']:
            pytest.xfail('not supported')
        if reference_model:
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


@pytest.fixture(scope="module", params=trials, ids=trial_names)
def io_trial(request, data_directory):
    reference_model = None
    if request.param.reference_file:
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



