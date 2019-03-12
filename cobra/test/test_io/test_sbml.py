# -*- coding: utf-8 -*-
"""
Testing SBML functionality based on libsbml.
"""

from __future__ import absolute_import

from collections import namedtuple
from os import unlink
from os.path import join, split
from pickle import load
from tempfile import gettempdir

import pytest

from cobra import Model
from cobra.io import read_sbml_model, validate_sbml_model, write_sbml_model


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
                  read_sbml_model, write_sbml_model,
                  validate_sbml_model),
          IOTrial('fbc2Gz', 'mini.pickle', 'mini_fbc2.xml.gz',
                  read_sbml_model, write_sbml_model, None),
          IOTrial('fbc2Bz2', 'mini.pickle', 'mini_fbc2.xml.bz2',
                  read_sbml_model, write_sbml_model, None),
          IOTrial('fbc1', 'mini.pickle', 'mini_fbc1.xml',
                  read_sbml_model, write_sbml_model, None),
          IOTrial('cobra', None, 'mini_cobra.xml',
                  read_sbml_model, write_sbml_model, None),
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
        assert abs(solution1.objective_value -
                   solution2.objective_value) < 0.001
        # ensure the references are correct
        assert model2.metabolites[0]._model is model2
        assert model2.reactions[0]._model is model2
        assert model2.genes[0]._model is model2

    @classmethod
    def extra_comparisons(cls, name, model1, model2):
        assert model1.compartments == model2.compartments

        # FIXME: problems of duplicate annotations in test data
        #  ('cas': ['56-65-5', '56-65-5'])
        # assert dict(model1.metabolites[4].annotation) == dict(
        #    model2.metabolites[4].annotation)
        d1 = model1.reactions[4].annotation
        d2 = model2.reactions[4].annotation
        assert list(d1.keys()) == list(d2.keys())
        for k in d1:
            assert set(d1[k]) == set(d2[k])
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


def test_filehandle(data_directory, tmp_path):
    """Test reading and writing to file handle."""
    with open(join(data_directory, "mini_fbc2.xml"), "r") as f_in:
        model1 = read_sbml_model(f_in)
        assert model1 is not None

    sbml_path = join(str(tmp_path), "test.xml")
    with open(sbml_path, "w") as f_out:
        write_sbml_model(model1, f_out)

    with open(sbml_path, "r") as f_in:
        model2 = read_sbml_model(f_in)

    TestCobraIO.compare_models(name="filehandle",
                               model1=model1, model2=model2)


def test_from_sbml_string(data_directory):
    """Test reading from SBML string."""
    sbml_path = join(data_directory, "mini_fbc2.xml")
    with open(sbml_path, "r") as f_in:
        sbml_str = f_in.read()
        model1 = read_sbml_model(sbml_str)

    model2 = read_sbml_model(sbml_path)
    TestCobraIO.compare_models(name="read from string",
                               model1=model1, model2=model2)


def test_model_history(tmp_path):
    """Testing reading and writing of ModelHistory."""
    model = Model("test")
    model._sbml = {
        "creators": [{
            "familyName": "Mustermann",
            "givenName": "Max",
            "organisation": "Muster University",
            "email": "muster@university.com",
        }]
    }

    sbml_path = join(str(tmp_path), "test.xml")
    with open(sbml_path, "w") as f_out:
        write_sbml_model(model, f_out)

    with open(sbml_path, "r") as f_in:
        model2 = read_sbml_model(f_in)

    assert "creators" in model2._sbml
    assert len(model2._sbml["creators"]) is 1
    c = model2._sbml["creators"][0]
    assert c["familyName"] == "Mustermann"
    assert c["givenName"] == "Max"
    assert c["organisation"] == "Muster University"
    assert c["email"] == "muster@university.com"


def test_groups(data_directory, tmp_path):
    """Testing reading and writing of groups"""
    sbml_path = join(data_directory, "e_coli_core.xml")
    model = read_sbml_model(sbml_path)
    assert model.groups is not None
    assert len(model.groups) == 10
    g1 = model.groups[0]
    assert len(g1.members) == 6

    temp_path = join(str(tmp_path), "test.xml")
    with open(temp_path, "w") as f_out:
        write_sbml_model(model, f_out)

    with open(temp_path, "r") as f_in:
        model2 = read_sbml_model(f_in)

        assert model2.groups is not None
        assert len(model2.groups) == 10
        g1 = model2.groups[0]
        assert len(g1.members) == 6


def test_validate(data_directory):
    """Test the validation code. """
    sbml_path = join(data_directory, "mini_fbc2.xml")
    with open(sbml_path, "r") as f_in:
        model1, errors = validate_sbml_model(f_in,
                                             check_modeling_practice=True)
        assert model1
        assert errors
        assert len(errors["SBML_WARNING"]) == 23


def test_infinity_bounds(data_directory, tmp_path):
    """Test infinity bound example. """
    sbml_path = join(data_directory, "fbc_ex1.xml")
    model = read_sbml_model(sbml_path)

    # check that simulation works
    solution = model.optimize()

    # check that values are set
    r = model.reactions.get_by_id("EX_X")
    assert r.lower_bound == -float("Inf")
    assert r.upper_bound == float("Inf")

    temp_path = join(str(tmp_path), "test.xml")
    with open(temp_path, "w") as f_out:
        write_sbml_model(model, f_out)

    with open(temp_path, "r") as f_in:
        model2 = read_sbml_model(f_in)
        r = model2.reactions.get_by_id("EX_X")
        assert r.lower_bound == -float("Inf")
        assert r.upper_bound == float("Inf")


def test_boundary_conditions(data_directory):
    """Test infinity bound example. """
    sbml_path1 = join(data_directory, "fbc_ex1.xml")
    model1 = read_sbml_model(sbml_path1)
    sol1 = model1.optimize()

    # model with species boundaryCondition==True
    sbml_path2 = join(data_directory, "fbc_ex2.xml")
    model2 = read_sbml_model(sbml_path2)
    sol2 = model2.optimize()

    r = model2.reactions.get_by_id("EX_X")
    assert r.lower_bound == -float("Inf")
    assert r.upper_bound == float("Inf")

    assert sol1.objective_value == sol2.objective_value
