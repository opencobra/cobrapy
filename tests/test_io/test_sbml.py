"""Testing SBML functionality based on libsbml."""

from collections import namedtuple
from os import unlink
from os.path import join, split
from pathlib import Path
from pickle import load
from tempfile import gettempdir
from typing import List, Tuple

import pytest
from _pytest.fixtures import SubRequest

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
IOTrial = namedtuple(
    "IOTrial",
    [
        "name",
        "reference_file",
        "test_file",
        "read_function",
        "write_function",
        "validation_function",
    ],
)
trials: List[IOTrial] = [
    IOTrial(
        "fbc2",
        "mini.pickle",
        "mini_fbc2.xml",
        read_sbml_model,
        write_sbml_model,
        validate_sbml_model,
    ),
    IOTrial(
        "fbc2Gz",
        "mini.pickle",
        "mini_fbc2.xml.gz",
        read_sbml_model,
        write_sbml_model,
        None,
    ),
    IOTrial(
        "fbc2Bz2",
        "mini.pickle",
        "mini_fbc2.xml.bz2",
        read_sbml_model,
        write_sbml_model,
        None,
    ),
    IOTrial(
        "fbc1", "mini.pickle", "mini_fbc1.xml", read_sbml_model, write_sbml_model, None
    ),
    IOTrial("cobra", None, "mini_cobra.xml", read_sbml_model, write_sbml_model, None),
]
trial_names: list = [node.name for node in trials]


@pytest.mark.parametrize("trial", trials)
def test_validate(trial: IOTrial, data_directory: Path) -> None:
    """Test validation function.

    Parameters
    ----------
    trial: IOTrial
        Which model trial to check.
    data_directory: Path
        Directory where the data is.
    """
    if trial.validation_function is None:
        pytest.skip("not implemented")
    test_file = data_directory / trial.test_file
    trial.validation_function(test_file)


class TestCobraIO:
    """Tests the read and write functions."""

    @classmethod
    def compare_models(cls, name: str, model1: Model, model2: Model) -> None:
        """Compare two models.

        name, str
            name of models to compare.
        model1: Model
            First model to compare.
        model2: Model
            Second model to compare.
        """
        print(name)
        assert len(model1.reactions) == len(model2.reactions)
        assert len(model1.metabolites) == len(model2.metabolites)
        assert model1.objective.direction == model2.objective.direction
        for attr in (
            "id",
            "name",
            "lower_bound",
            "upper_bound",
            "objective_coefficient",
            "gene_reaction_rule",
        ):
            assert getattr(model1.reactions[0], attr) == getattr(
                model2.reactions[0], attr
            )
            assert getattr(model1.reactions[5], attr) == getattr(
                model2.reactions[5], attr
            )
            assert getattr(model1.reactions[-1], attr) == getattr(
                model2.reactions[-1], attr
            )
        for attr in ("id", "name", "compartment", "formula", "charge"):
            assert getattr(model1.metabolites[0], attr) == getattr(
                model2.metabolites[0], attr
            )
            assert getattr(model1.metabolites[5], attr) == getattr(
                model2.metabolites[5], attr
            )
            assert getattr(model1.metabolites[-1], attr) == getattr(
                model2.metabolites[-1], attr
            )
        assert len(model1.reactions[0].metabolites) == len(
            model2.reactions[0].metabolites
        )
        assert len(model1.reactions[8].metabolites) == len(
            model2.reactions[8].metabolites
        )
        assert len(model1.reactions[-1].metabolites) == len(
            model2.reactions[-1].metabolites
        )
        assert len(model1.genes) == len(model2.genes)

        # ensure they have the same solution max
        solution1 = model1.optimize()
        solution2 = model2.optimize()
        assert (
            solution1.status == "infeasible" and solution2.status == "infeasible"
        ) or abs(solution1.objective_value - solution2.objective_value) < 0.001
        # ensure the references are correct
        assert model2.metabolites[0]._model is model2
        assert model2.reactions[0]._model is model2
        assert model2.genes[0]._model is model2

    @classmethod
    def extra_comparisons(cls, name: str, model1: Model, model2: Model) -> None:
        """Compare additional features of the model.

        name, str
            name of models to compare.
        model1: Model
            First model to compare.
        model2: Model
            Second model to compare.
        """
        print(name)
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
            model2.reactions[4].annotation
        )
        assert dict(model1.genes[5].annotation) == dict(model2.genes[5].annotation)

        for attr in ("id", "name"):
            assert getattr(model1.genes[0], attr) == getattr(model2.genes[0], attr)
            assert getattr(model1.genes[10], attr) == getattr(model2.genes[10], attr)
            assert getattr(model1.genes[-1], attr) == getattr(model2.genes[-1], attr)

    def test_read_1(self, io_trial: Tuple[str, Model, Model, Model]) -> None:
        """Read the first model from a processed IOTrial via io_trial().

        io_trial: [str, Model, Model, Model]
            A tuple containing name, reference_model, test_model, rewritten_model.
        """
        name, reference_model, test_model, _ = io_trial
        if name in ["fbc1"]:
            pytest.xfail("not supported")
        if reference_model:
            self.compare_models(name, reference_model, test_model)

    def test_read_2(self, io_trial: Tuple[str, Model, Model, Model]) -> None:
        """Read the second model from a processed IOTrial via io_trial().

        io_trial: [str, Model, Model, Model]
            A tuple containing name, reference_model, test_model, rewritten_model.
        """
        name, reference_model, test_model, _ = io_trial
        if name in ["fbc1", "mat", "cobra", "raven-mat"]:
            pytest.xfail("not supported")
        if reference_model:
            self.extra_comparisons(name, reference_model, test_model)

    def test_write_1(self, io_trial: Tuple[str, Model, Model, Model]) -> None:
        """Test writing the first model from a processed IOTrial via io_trial().

        io_trial: [str, Model, Model, Model]
            A tuple containing name, reference_model, test_model, rewritten_model.
        """
        name, _, test_model, reread_model = io_trial
        if name in ["fbc1", "raven-mat"]:
            pytest.xfail("not supported")

        self.compare_models(name, test_model, reread_model)

    def test_write_2(self, io_trial: Tuple[str, Model, Model, Model]) -> None:
        """Test writing the second model from a processed IOTrial via io_trial().

        io_trial: [str, Model, Model, Model]
            A tuple containing name, reference_model, test_model, rewritten_model.
        """
        name, _, test_model, reread_model = io_trial
        if name in ["fbc1", "mat", "cobra", "raven-mat"]:
            pytest.xfail("not supported")
        self.extra_comparisons(name, test_model, reread_model)


@pytest.fixture(scope="module", params=trials, ids=trial_names)
def io_trial(
    request: SubRequest, data_directory: Path
) -> Tuple[str, Model, Model, Model]:
    """Read reference model, test model, write test model and reread it.

    Parameters
    ----------
    request: IOTrail
    data_directory: Path
        Directory where the data is.

    This function will read the reference model, the test model. It will then write
    the test model based on the write_function() in IOTrial, and read the written test
    model using read_function in IOTrail.
    This can be used to compare the original model, test model and the written and
    reread model.

    Returns
    -------
    Tuple: str, Model, Model, Model
        Name, original model as read by read_function(), test model as read by
        read_function(), test model as written by write_function() and reread by
        read_function().
    """
    reference_model = None
    if request.param.reference_file:
        with open(
            data_directory.joinpath(request.param.reference_file), "rb"
        ) as infile:
            reference_model = load(infile)
    test_model = request.param.read_function(data_directory / request.param.test_file)
    test_output_filename = join(gettempdir(), split(request.param.test_file)[-1])
    # test writing the model within a context with a non-empty stack
    with test_model:
        test_model.objective = test_model.objective
        request.param.write_function(test_model, test_output_filename)
    reread_model = request.param.read_function(test_output_filename)
    unlink(test_output_filename)
    return request.param.name, reference_model, test_model, reread_model


def test_filehandle(data_directory: Path, tmp_path: Path) -> None:
    """Test reading and writing to file handle.

    Parameters
    ----------
    data_directory: Path
        Directory where the data is.
    tmp_path: Path
        Directory to use for temporary data.
    """
    with data_directory.joinpath("mini_fbc2.xml").open("r") as f_in:
        model1 = read_sbml_model(f_in)
        assert model1 is not None

    sbml_path = tmp_path / "test.xml"
    with sbml_path.open("w") as f_out:
        write_sbml_model(model1, f_out)

    with sbml_path.open("r") as f_in:
        model2 = read_sbml_model(f_in)

    TestCobraIO.compare_models(name="filehandle", model1=model1, model2=model2)


def test_from_sbml_string(data_directory: Path) -> None:
    """Test reading from SBML string.

    Parameters
    ----------
    data_directory: Path
        Directory where the data is.
    """
    sbml_path = data_directory / "mini_fbc2.xml"
    model1 = read_sbml_model(sbml_path.read_text())

    model2 = read_sbml_model(sbml_path)
    TestCobraIO.compare_models(name="read from string", model1=model1, model2=model2)


@pytest.mark.skip(reason="Model history currently not written")
def test_model_history(tmp_path: Path) -> None:
    """Testing reading and writing of ModelHistory.

    Parameters
    ----------
    tmp_path: Path
        Directory to use for temporary data.
    """
    model = Model("test")
    model._sbml = {
        "creators": [
            {
                "familyName": "Mustermann",
                "givenName": "Max",
                "organisation": "Muster University",
                "email": "muster@university.com",
            }
        ]
    }

    sbml_path = tmp_path / "test.xml"
    with open(sbml_path, "w") as f_out:
        write_sbml_model(model, f_out)

    with open(sbml_path, "r") as f_in:
        model2 = read_sbml_model(f_in)

    assert "creators" in model2._sbml
    assert len(model2._sbml["creators"]) == 1
    c = model2._sbml["creators"][0]
    assert c["familyName"] == "Mustermann"
    assert c["givenName"] == "Max"
    assert c["organisation"] == "Muster University"
    assert c["email"] == "muster@university.com"


def test_groups(data_directory: Path, tmp_path: Path) -> None:
    """Testing reading and writing of groups.

    Parameters
    ----------
    data_directory: Path
        Directory where the data is.
    tmp_path: Path
        Directory to use for temporary data.
    """
    sbml_path = data_directory / "e_coli_core.xml"
    model = read_sbml_model(sbml_path)
    assert model.groups is not None
    assert len(model.groups) == 10
    g1 = model.groups[0]
    assert len(g1.members) == 6

    temp_path = tmp_path / "test.xml"
    with open(temp_path, "w") as f_out:
        write_sbml_model(model, f_out)

    with open(temp_path, "r") as f_in:
        model2 = read_sbml_model(f_in)

        assert model2.groups is not None
        assert len(model2.groups) == 10
        g1 = model2.groups[0]
        assert len(g1.members) == 6


def test_missing_flux_bounds1(data_directory: Path) -> None:
    """Test missing flux bounds in an incorrect model.

    Parameters
    ----------
    data_directory: Path
        Directory where the data is.
    """
    sbml_path = data_directory / "annotation.xml"
    with open(sbml_path, "r") as f_in:
        # missing flux bounds are set to cobra.configuration.bounds
        # noinspection PyTupleAssignmentBalance
        model, errors = validate_sbml_model(f_in, set_missing_bounds=True)
        r1 = model.reactions.R1
        assert r1.lower_bound == config.lower_bound
        assert r1.upper_bound == config.upper_bound


def test_missing_flux_bounds2(data_directory: Path) -> None:
    """Test missing flux bounds set to [-INF, INF].

    Parameters
    ----------
    data_directory: Path
        Directory where the data is.
    """
    sbml_path = data_directory / "annotation.xml"
    with open(sbml_path, "r") as f_in:
        # missing flux bounds are set to [-INF, INF]
        # noinspection PyTupleAssignmentBalance
        model, errors = validate_sbml_model(f_in, set_missing_bounds=False)
        r1 = model.reactions.R1
        assert r1.lower_bound == config.lower_bound
        assert r1.upper_bound == config.upper_bound


def test_validate2(data_directory: Path) -> None:
    """Test the validation code.

    Parameters
    ----------
    data_directory: Path
        Directory where the data is.
    """
    sbml_path = data_directory / "mini_fbc2.xml"
    with open(sbml_path, "r") as f_in:
        # noinspection PyTupleAssignmentBalance
        model1, errors = validate_sbml_model(f_in, check_modeling_practice=True)
        assert model1
        assert errors
        assert len(errors["SBML_WARNING"]) == 0


def test_validation_warnings(data_directory: Path) -> None:
    """Test the validation warnings.

    Parameters
    ----------
    data_directory: Path
        Directory where the data is.
    """
    sbml_path = data_directory / "validation.xml"
    with open(sbml_path, "r") as f_in:
        # noinspection PyTupleAssignmentBalance
        model1, errors = validate_sbml_model(f_in, check_modeling_practice=True)
        assert model1
        assert errors
        assert len(errors["COBRA_WARNING"]) == 3
        assert "No objective in listOfObjectives" in errors["COBRA_WARNING"]


def test_infinity_bounds(data_directory: Path, tmp_path: Path) -> None:
    """Test infinity bound example.

    Parameters
    ----------
    data_directory: Path
        Directory where the data is.
    tmp_path: Path
        Directory to use for temporary data.
    """
    sbml_path = data_directory / "fbc_ex1.xml"
    model = read_sbml_model(sbml_path)

    # check that simulation works
    solution = model.optimize()
    assert solution is not None

    # check that values are set
    r = model.reactions.get_by_id("EX_X")
    assert r.lower_bound == -float("Inf")
    assert r.upper_bound == float("Inf")

    temp_path = tmp_path / "test.xml"
    with open(temp_path, "w") as f_out:
        write_sbml_model(model, f_out)

    with open(temp_path, "r") as f_in:
        model2 = read_sbml_model(f_in)
        r = model2.reactions.get_by_id("EX_X")
        assert r.lower_bound == -float("Inf")
        assert r.upper_bound == float("Inf")


def test_boundary_conditions(data_directory: Path) -> None:
    """Test infinity bound example.

    Parameters
    ----------
    data_directory: Path
        Directory where the data is.
    """
    sbml_path1 = data_directory / "fbc_ex1.xml"
    model1 = read_sbml_model(sbml_path1)
    sol1 = model1.optimize()

    # model with species boundaryCondition==True
    sbml_path2 = data_directory / "fbc_ex2.xml"
    model2 = read_sbml_model(sbml_path2)
    sol2 = model2.optimize()

    r = model2.reactions.get_by_id("EX_X")
    assert r.lower_bound == config.lower_bound
    assert r.upper_bound == config.upper_bound

    assert sol1.objective_value == sol2.objective_value


def test_bounds_on_write(data_directory: Path, tmp_path: Path) -> None:
    """Test infinity bound example.

    Parameters
    ----------
    data_directory: Path
        Directory where the data is.
    """
    sbml_path1 = data_directory / "fbc_ex1.xml"
    model1 = read_sbml_model(sbml_path1)

    r_x = model1.reactions.get_by_id("EX_X")
    r_y = model1.reactions.get_by_id("EX_Ac")

    r_x.bounds = (config.lower_bound - 1000, config.upper_bound + 1000)
    assert r_x.lower_bound == config.lower_bound - 1000
    assert r_x.upper_bound == config.upper_bound + 1000

    # Global min/max bounds for other reactions should not change before & after write!
    r_y.bounds = (config.lower_bound, config.upper_bound)
    assert r_y.lower_bound == config.lower_bound
    assert r_y.upper_bound == config.upper_bound

    sbml_path = tmp_path / "test.xml"
    with open(sbml_path, "w") as f_out:
        write_sbml_model(model1, f_out)

    with open(sbml_path, "r") as f_in:
        model2 = read_sbml_model(f_in)

    r2_x = model2.reactions.get_by_id("EX_X")
    r2_y = model2.reactions.get_by_id("EX_Ac")

    assert r2_x.lower_bound == r_x.lower_bound
    assert r2_x.upper_bound == r_x.upper_bound
    assert r2_y.lower_bound == r_y.lower_bound  # before fix #1300, this would fail
    assert r2_y.upper_bound == r_y.upper_bound  # before fix #1300, this would fail


def test_gprs(large_model: Model, tmp_path: Path) -> None:
    """Test that GPRs are written and read correctly.

    Parameters
    ----------
    large_model: Model
        Model to test gprs on.
    tmp_path: Path
        Directory to use for temporary data.
    """
    model1 = large_model
    sbml_path = tmp_path / "test.xml"
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


def test_identifiers_annotation() -> None:
    """Test annotation with identifiers."""
    from cobra.io.sbml import _parse_annotation_info

    for uri in [
        "http://identifiers.org/chebi/CHEBI:000123",
        "https://identifiers.org/chebi/CHEBI:000123",
        "http://identifiers.org/CHEBI:000123",
        "https://identifiers.org/CHEBI:000123",
    ]:
        data = _parse_annotation_info(uri)
        assert data
        assert data[0] == "chebi"
        assert data[1] == "CHEBI:000123"

    for uri in [
        "http://identifiers.org/taxonomy/9602",
        "https://identifiers.org/taxonomy/9602",
        "http://identifiers.org/taxonomy:9602",
        "https://identifiers.org/taxonomy:9602",
    ]:
        data = _parse_annotation_info(uri)
        assert data
        assert data[0] == "taxonomy"
        assert data[1] == "9602"

    for uri in [
        "http://identifier.org/taxonomy/9602",
        "https://test.com",
    ]:
        data = _parse_annotation_info(uri)
        assert data is None


def test_smbl_with_notes(data_directory: Path, tmp_path: Path) -> None:
    """Test that NOTES in the RECON 2.2 style are written and read correctly.

    Parameters
    ----------
    data_directory: Path
        Directory where the data is.
    tmp_path: Path
        Directory to use for temporary data.
    """
    sbml_path = join(data_directory, "example_notes.xml")
    model = read_sbml_model(sbml_path)
    assert model.metabolites is not None
    metabolite_notes = {
        "2hb_e": {"CHARGE": "-1", "FORMULA": "C4H7O3", "SMILES": "CCC(O)C(O)=O"},
        "nad_e": {
            "CHARGE": "-1",
            "FORMULA": "C21H26N7O14P2",
            "SMILES": "NC(=O)c1ccc[n+](c1)[C@@H]1O[C@H](COP([O-])("
            "=O)OP([O-])(=O)OC[C@H]2O[C@H]([C@H](O)["
            "C@@H]2O)n2cnc3c(N)ncnc23)[C@@H](O)[C@H]1O",
        },
        "h_e": {"CHARGE": "1", "FORMULA": "H", "SMILES": "[1H+]"},
        "2obut_e": {"CHARGE": "-1", "FORMULA": "C4H5O3", "SMILES": "CCC(=O)C([O-])=O"},
        "nadh_e": {
            "CHARGE": "-2",
            "FORMULA": "C21H27N7O14P2",
            "SMILES": "NC(=O)C1=CN(C=CC1)[C@@H]1O[C@H](COP([O-])("
            "=O)OP([O-])(=O)OC[C@H]2O[C@H]([C@H](O)["
            "C@@H]2O)n2cnc3c(N)ncnc23)[C@@H](O)[C@H]1O",
        },
    }
    metabolite_annotations = {
        "2hb_e": {
            "sbo": "SBO:0000247",
            "inchi": "InChI=1S/C4H8O3/c1-2-3(5)4(6)7/h3,5H,2H2,1H3,(H,6,7)",
            "chebi": "CHEBI:1148",
        },
        "nad_e": {
            "sbo": "SBO:0000247",
            "inchi": "InChI=1S/C21H27N7O14P2/c22-17-12-19("
            "25-7-24-17)28(8-26-12)21-16(32)14(30)11("
            "41-21)6-39-44(36,37)42-43(34,35)38-5-10-13(29)15("
            "31)20(40-10)27-3-1-2-9(4-27)18("
            "23)33/h1-4,7-8,10-11,13-16,20-21,29-32H,5-6H2,"
            "(H5-,22,23,24,25,33,34,35,36,37)/p-1/t10-,"
            "11-,13-,14-,15-,16-,20-,21-/m1/s1",
            "chebi": "CHEBI:57540",
        },
        "h_e": {
            "sbo": "SBO:0000247",
            "inchi": "InChI=1S/p+1/i/hH",
            "chebi": "CHEBI:24636",
        },
        "2obut_e": {
            "sbo": "SBO:0000247",
            "inchi": "InChI=1S/C4H6O3/c1-2-3(5)4(6)7/h2H2,1H3,(H,6,7)/p-1",
            "chebi": "CHEBI:16763",
        },
        "nadh_e": {
            "sbo": "SBO:0000247",
            "inchi": "InChI=1S/C21H29N7O14P2/c22-17-12-19("
            "25-7-24-17)28(8-26-12)21-16(32)14(30)11("
            "41-21)6-39-44(36,37)42-43(34,35)38-5-10-13("
            "29)15(31)20(40-10)27-3-1-2-9(4-27)18("
            "23)33/h1,3-4,7-8,10-11,13-16,20-21,29-32H,2,"
            "5-6H2,(H2,23,33)(H,34,35)(H,36,37)(H2,22,24,"
            "25)/p-2/t10-,11-,13-,14-,15-,16-,20-,21-/m1/s1",
            "chebi": "CHEBI:57945",
        },
    }
    reaction_notes = {
        "CONFIDENCE_LEVEL": "4",
        "NOTES": "NCD",
        "SUBSYSTEM": "Propanoate metabolism",
        "GENE_ASSOCIATION": "(HGNC:8546 and HGNC:8548) or (HGNC:8547 and HGNC:8548)",
    }
    reaction_annotations = {
        "sbo": "SBO:0000176",
        "ec-code": "1.1.1.27",
        "pubmed": ["10108", "21765"],
    }

    for met_id in metabolite_notes:
        assert model.metabolites.has_id(met_id)
        for note_key in metabolite_notes[met_id].keys():
            assert note_key in model.metabolites.get_by_id(met_id).notes
            assert (
                metabolite_notes[met_id][note_key]
                == model.metabolites.get_by_id(met_id).notes[note_key]
            )
        for annotation_key in metabolite_annotations[met_id].keys():
            assert annotation_key in model.metabolites.get_by_id(met_id).annotation
            assert (
                metabolite_annotations[met_id][annotation_key]
                == model.metabolites.get_by_id(met_id).annotation[annotation_key]
            )

    for note_key in reaction_notes.keys():
        assert note_key in model.reactions[0].notes.keys()
        assert reaction_notes[note_key] == model.reactions[0].notes[note_key]
    assert (
        model.reactions[0].gene_reaction_rule == "(HGNC:8546 and HGNC:8548) or "
        "(HGNC:8547 and HGNC:8548)"
    )
    assert len(model.groups) == 1
    for annotation_key in reaction_annotations.keys():
        assert annotation_key in model.reactions[0].annotation.keys()
        assert (
            reaction_annotations[annotation_key]
            == model.reactions[0].annotation[annotation_key]
        )


def test_stable_gprs(data_directory: Path, tmp_path: Path) -> None:
    """Test that GPRs are written correctly after manual changes.

    Parameters
    ----------
    data_directory: Path
        Directory where the data is.
    tmp_path: Path
        Directory to use for temporary data.
    """
    mini = read_sbml_model(join(data_directory, "mini_fbc2.xml"))
    mini.reactions.GLCpts.gene_reaction_rule = "((b2415 and b2417)or (b2416))"
    fixed = join(str(tmp_path), "fixed_gpr.xml")
    write_sbml_model(mini, fixed)
    fixed_model = read_sbml_model(fixed)
    assert (
        fixed_model.reactions.GLCpts.gene_reaction_rule == "(b2415 and b2417) or b2416"
    )


def test_history(data_directory: Path) -> None:
    """Test that the history is read from the model."""
    mini = read_sbml_model(join(data_directory, "mini_history.xml"))
    assert "creators" in mini._sbml
    assert "organisation" in mini._sbml["creators"][0]
    assert "created" in mini._sbml
    assert isinstance(mini._sbml["created"], str)
    assert "2022" in mini._sbml["created"]
