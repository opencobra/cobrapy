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

import cobra
from cobra import Model
from cobra.io import read_sbml_model, validate_sbml_model, write_sbml_model


config = cobra.Configuration()  # for default bounds

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
        assert (solution1.status == 'infeasible' and
                solution2.status == 'infeasible') or \
            abs(solution1.objective_value -
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


@pytest.mark.skip(reason="Model history currently not written")
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


def test_missing_flux_bounds1(data_directory):
    sbml_path = join(data_directory, "annotation.xml")
    with open(sbml_path, "r") as f_in:
        # missing flux bounds are set to cobra.configuration.bounds
        model, errors = validate_sbml_model(f_in,
                                            set_missing_bounds=True)
        r1 = model.reactions.R1
        assert r1.lower_bound == config.lower_bound
        assert r1.upper_bound == config.upper_bound


def test_missing_flux_bounds2(data_directory):
    sbml_path = join(data_directory, "annotation.xml")
    with open(sbml_path, "r") as f_in:
        # missing flux bounds are set to [-INF, INF]
        model, errors = validate_sbml_model(f_in,
                                            set_missing_bounds=False)
        r1 = model.reactions.R1
        assert r1.lower_bound == config.lower_bound
        assert r1.upper_bound == config.upper_bound


def test_validate(data_directory):
    """Test the validation code. """
    sbml_path = join(data_directory, "mini_fbc2.xml")
    with open(sbml_path, "r") as f_in:
        model1, errors = validate_sbml_model(f_in,
                                             check_modeling_practice=True)
        assert model1
        assert errors
        assert len(errors["SBML_WARNING"]) == 0


def test_validation_warnings(data_directory):
    """Test the validation warnings. """
    sbml_path = join(data_directory, "validation.xml")
    with open(sbml_path, "r") as f_in:
        model1, errors = validate_sbml_model(f_in,
                                             check_modeling_practice=True)
        assert model1
        assert errors
        assert len(errors["COBRA_WARNING"]) == 3
        assert "No objective in listOfObjectives" in errors["COBRA_WARNING"]


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
    assert r.lower_bound == config.lower_bound
    assert r.upper_bound == config.upper_bound

    assert sol1.objective_value == sol2.objective_value


def test_gprs(data_directory, tmp_path):
    """Test that GPRs are written and read correctly"""
    model1 = read_sbml_model(join(data_directory, "iJO1366.xml.gz"))

    sbml_path = join(str(tmp_path), "test.xml")
    with open(sbml_path, "w") as f_out:
        write_sbml_model(model1, f_out)

    with open(sbml_path, "r") as f_in:
        model2 = read_sbml_model(f_in)

    for r1 in model1.reactions:
        rid = r1.id
        r2 = model2.reactions.get_by_id(rid)
        gpr1 = r1.gene_reaction_rule
        gpr2 = r2.gene_reaction_rule

        assert gpr1 == gpr2


def test_identifiers_annotation():
    from cobra.io.sbml import parse_annotation_info

    for uri in [
        "http://identifiers.org/chebi/CHEBI:000123",
        "https://identifiers.org/chebi/CHEBI:000123",
        "http://identifiers.org/CHEBI:000123",
        "https://identifiers.org/CHEBI:000123",
    ]:
        data = parse_annotation_info(uri)
        assert data
        assert data[0] == "chebi"
        assert data[1] == "CHEBI:000123"

    for uri in [
        "http://identifiers.org/taxonomy/9602",
        "https://identifiers.org/taxonomy/9602",
        "http://identifiers.org/taxonomy:9602",
        "https://identifiers.org/taxonomy:9602",
    ]:
        data = parse_annotation_info(uri)
        assert data
        assert data[0] == "taxonomy"
        assert data[1] == "9602"

    for uri in [
        "http://identifier.org/taxonomy/9602",
        "https://test.com",
    ]:
        data = parse_annotation_info(uri)
        assert data is None


def test_smbl_with_notes(data_directory, tmp_path):
    """Test that NOTES in the RECON 2.2 style are written and read correctly"""
    sbml_path = join(data_directory, "example_notes.xml")
    model = read_sbml_model(sbml_path)
    assert model.metabolites is not None
    metabolite_notes = {
        '2hb_e': {'CHARGE': '-1', 'FORMULA': 'C4H7O3',
                  'SMILES': 'CCC(O)C(O)=O'},
        'nad_e': {'CHARGE': '-1', 'FORMULA': 'C21H26N7O14P2',
                  'SMILES': 'NC(=O)c1ccc[n+](c1)[C@@H]1O[C@H](COP([O-])('
                            '=O)OP([O-])(=O)OC[C@H]2O[C@H]([C@H](O)['
                            'C@@H]2O)n2cnc3c(N)ncnc23)[C@@H](O)[C@H]1O'},
        'h_e': {'CHARGE': '1', 'FORMULA': 'H', 'SMILES': '[1H+]'},
        '2obut_e': {'CHARGE': '-1', 'FORMULA': 'C4H5O3',
                    'SMILES': 'CCC(=O)C([O-])=O'},
        'nadh_e': {'CHARGE': '-2', 'FORMULA': 'C21H27N7O14P2',
                   'SMILES': 'NC(=O)C1=CN(C=CC1)[C@@H]1O[C@H](COP([O-])('
                             '=O)OP([O-])(=O)OC[C@H]2O[C@H]([C@H](O)['
                             'C@@H]2O)n2cnc3c(N)ncnc23)[C@@H](O)[C@H]1O'}
    }
    metabolite_annotations = {
        '2hb_e': {'sbo': 'SBO:0000247',
                  'inchi': ['InChI=1S/C4H8O3/c1-2-3(5)4(6)7/h3,5H,2H2,1H3,'
                            '(H,6,7)'],
                  'chebi': ['CHEBI:1148']},
        'nad_e': {'sbo': 'SBO:0000247',
                  'inchi': ['InChI=1S/C21H27N7O14P2/c22-17-12-19('
                            '25-7-24-17)28(8-26-12)21-16(32)14(30)11('
                            '41-21)6-39-44(36,37)42-43(34,35)38-5-10-13(29)15('
                            '31)20(40-10)27-3-1-2-9(4-27)18('
                            '23)33/h1-4,7-8,10-11,13-16,20-21,29-32H,5-6H2,'
                            '(H5-,22,23,24,25,33,34,35,36,37)/p-1/t10-,'
                            '11-,13-,14-,15-,16-,20-,21-/m1/s1'],
                  'chebi': ['CHEBI:57540']},
        'h_e': {'sbo': 'SBO:0000247', 'inchi': ['InChI=1S/p+1/i/hH'],
                'chebi': ['CHEBI:24636']},
        '2obut_e': {'sbo': 'SBO:0000247',
                    'inchi': ['InChI=1S/C4H6O3/c1-2-3(5)4(6)7/h2H2,1H3,(H,6,'
                              '7)/p-1'],
                    'chebi': ['CHEBI:16763']},
        'nadh_e': {'sbo': 'SBO:0000247',
                   'inchi': ['InChI=1S/C21H29N7O14P2/c22-17-12-19('
                             '25-7-24-17)28(8-26-12)21-16(32)14(30)11('
                             '41-21)6-39-44(36,37)42-43(34,35)38-5-10-13('
                             '29)15(31)20(40-10)27-3-1-2-9(4-27)18('
                             '23)33/h1,3-4,7-8,10-11,13-16,20-21,29-32H,2,'
                             '5-6H2,(H2,23,33)(H,34,35)(H,36,37)(H2,22,24,'
                             '25)/p-2/t10-,11-,13-,14-,15-,16-,20-,21-/m1/s1'],
                   'chebi': ['CHEBI:57945']}
    }
    reaction_notes = {'CONFIDENCE_LEVEL': '4', 'NOTES': 'NCD',
                      'SUBSYSTEM': 'Propanoate metabolism',
                      'GENE_ASSOCIATION': '(HGNC:8546 and HGNC:8548) or'
                                          ' (HGNC:8547 and HGNC:8548)'}
    reaction_annotations = {'sbo': 'SBO:0000176', 'ec-code': ['1.1.1.27'],
                            'pubmed': ['10108', '21765']}

    for met_id in metabolite_notes:
        assert model.metabolites.has_id(met_id)
        for note_key in metabolite_notes[met_id].keys():
            assert note_key in model.metabolites.get_by_id(met_id).notes
            assert metabolite_notes[met_id][note_key] == \
                model.metabolites.get_by_id(met_id).notes[note_key]
        for annotation_key in metabolite_annotations[met_id].keys():
            assert annotation_key in model.metabolites.get_by_id(
                met_id).annotation
            print(met_id)
            assert metabolite_annotations[met_id][annotation_key] == \
                model.metabolites.get_by_id(
                    met_id).annotation[annotation_key]

    for note_key in reaction_notes.keys():
        assert note_key in model.reactions[0].notes.keys()
        assert reaction_notes[note_key] == model.reactions[0].notes[note_key]
    assert model.reactions[
               0].gene_reaction_rule == '(HGNC:8546 and HGNC:8548) or ' \
                                        '(HGNC:8547 and HGNC:8548)'
    assert len(model.groups) == 1
    for annotation_key in reaction_annotations.keys():
        assert annotation_key in model.reactions[0].annotation.keys()
        assert reaction_annotations[annotation_key] == model.reactions[
            0].annotation[annotation_key]
