from warnings import catch_warnings, warn
from tempfile import gettempdir
from os import unlink
from os.path import join, split
from unittest import TestCase, TestLoader, TextTestRunner, skipIf, \
    expectedFailure
from functools import partial
from pickle import load, dump
import sys

if __name__ == "__main__":
    from cobra import io
    from cobra.test import data_directory
else:
    from .. import io
    from . import data_directory

libraries = ["scipy", "libsbml", "cPickle", "jsonschema"]
for library in libraries:
    try:
        exec("import %s" % library)
    except ImportError:
        exec("%s = None" % library)

with open(join(data_directory, "mini.pickle"), "rb") as infile:
    mini_model = load(infile)


class TestCobraIO(object):
    def compare_models(self, model1, model2):
        self.assertEqual(len(model1.reactions),
                         len(model2.reactions))
        self.assertEqual(len(model1.metabolites),
                         len(model2.metabolites))
        for attr in ("id", "name", "lower_bound", "upper_bound",
                     "objective_coefficient", "gene_reaction_rule"):
            self.assertEqual(getattr(model1.reactions[0], attr),
                             getattr(model2.reactions[0], attr))
            self.assertEqual(getattr(model1.reactions[10], attr),
                             getattr(model2.reactions[10], attr))
            self.assertEqual(getattr(model1.reactions[-1], attr),
                             getattr(model2.reactions[-1], attr))
        for attr in ("id", "name", "compartment", "formula", "charge"):
            self.assertEqual(getattr(model1.metabolites[0], attr),
                             getattr(model2.metabolites[0], attr))
            self.assertEqual(getattr(model1.metabolites[10], attr),
                             getattr(model2.metabolites[10], attr))
            self.assertEqual(getattr(model1.metabolites[-1], attr),
                             getattr(model2.metabolites[-1], attr))
        self.assertEqual(len(model1.reactions[0].metabolites),
                         len(model2.reactions[0].metabolites))
        self.assertEqual(len(model1.reactions[14].metabolites),
                         len(model2.reactions[14].metabolites))
        self.assertEqual(len(model1.reactions[-1].metabolites),
                         len(model2.reactions[-1].metabolites))
        self.assertEqual(len(model1.genes), len(model2.genes))
        # ensure they have the same solution max
        model1.optimize()
        model2.optimize()
        self.assertAlmostEqual(model1.solution.f, model2.solution.f,
                               places=3)
        # ensure the references are correct
        self.assertIs(model2.metabolites[0]._model, model2)
        self.assertIs(model2.reactions[0]._model, model2)
        self.assertIs(model2.genes[0]._model, model2)
        self.extra_comparisons(model1, model2)

    def extra_comparisons(self, model1, model2):
        # Overriding this prevents these extra comparisons. For example, mat
        # will need to skip this.
        self.assertEqual(model1.compartments, model2.compartments)
        self.assertEqual(model1.metabolites[4].annotation,
                         model2.metabolites[4].annotation)
        self.assertEqual(model1.reactions[4].annotation,
                         model2.reactions[4].annotation)
        self.assertEqual(model1.genes[4].annotation,
                         model2.genes[4].annotation)
        for attr in ("id", "name"):
            self.assertEqual(getattr(model1.genes[0], attr),
                             getattr(model2.genes[0], attr))
            self.assertEqual(getattr(model1.genes[10], attr),
                             getattr(model2.genes[10], attr))
            self.assertEqual(getattr(model1.genes[-1], attr),
                             getattr(model2.genes[-1], attr))
        return

    def test_read(self):
        read_model = self.read_function(self.test_file)
        self.compare_models(self.test_model, read_model)

    def test_read_nonexistent(self):
        # make sure that an error is raised when given a nonexistent file
        self.assertRaises(IOError, self.read_function, "fake_file")

    def test_write(self):
        test_output_filename = join(gettempdir(), split(self.test_file)[-1])
        self.write_function(self.test_model, test_output_filename)
        reread_model = self.read_function(test_output_filename)
        self.compare_models(self.test_model, reread_model)
        self.validate(test_output_filename)
        unlink(test_output_filename)

    def test_write_empty(self):
        test_output_filename = join(gettempdir(), split(self.test_file)[-1])
        m = self.test_model.copy()
        m.metabolites[0].charge = None
        m.remove_reactions(list(m.reactions))
        self.write_function(m, test_output_filename)
        reread_model = self.read_function(test_output_filename)
        self.assertEqual(len(reread_model.reactions), 0)
        self.assertEqual(len(reread_model.metabolites), len(m.metabolites))
        # ensure empty metabolite charge is read as None
        self.assertIs(reread_model.metabolites[0].charge, None)
        unlink(test_output_filename)

    def validate(self, filename):
        # overload if a validator exists
        None


class TestCobraIOSBMLfbc2(TestCase, TestCobraIO):
    def setUp(self):
        self.test_model = mini_model
        self.test_file = join(data_directory, "mini_fbc2.xml")
        self.read_function = io.read_sbml_model
        self.write_function = io.write_sbml_model

    def validate(self, filename):
        io.sbml3.validate_sbml_model(filename)


class TestCobraIOSBMLfbc2Gz(TestCobraIOSBMLfbc2):
    def setUp(self):
        self.test_model = mini_model
        self.test_file = join(data_directory, "mini_fbc2.xml.gz")
        self.read_function = io.read_sbml_model
        self.write_function = io.write_sbml_model


class TestCobraIOSBMLfbc2Bz2(TestCobraIOSBMLfbc2):
    def setUp(self):
        self.test_model = mini_model
        self.test_file = join(data_directory, "mini_fbc2.xml.bz2")
        self.read_function = io.read_sbml_model
        self.write_function = io.write_sbml_model


class TestCobraSBMLValidation(TestCase):
    def test_bad_valiation(self):
        for i in range(3):
            filename = join(data_directory, "invalid%d.xml" % i)
            m, errors = io.sbml3.validate_sbml_model(filename)
            self.assertTrue(len(errors) >= 1)


@skipIf(not libsbml, "libsbml required")
class TestCobraIOSBMLfbc1(TestCase, TestCobraIO):
    def setUp(self):
        self.test_model = mini_model
        self.test_file = join(data_directory, "mini_fbc1.xml")
        self.read_function = io.read_sbml_model
        self.write_function = partial(io.write_legacy_sbml,
                                      use_fbc_package=True)

    def extra_comparisons(self, model1, model2):
        None

    @expectedFailure
    def test_read(self):
        TestCobraIO.test_read(self)

    @expectedFailure
    def test_write(self):
        TestCobraIO.test_write(self)


@skipIf(not libsbml, "libsbml required")
class TestCobraIOSBMLcobra(TestCase, TestCobraIO):
    def setUp(self):
        self.test_model = mini_model
        self.test_file = join(data_directory, "mini_cobra.xml")
        self.read_function = io.read_sbml_model
        self.write_function = partial(io.write_sbml_model,
                                      use_fbc_package=False)

    def extra_comparisons(self, model1, model2):
        None


@skipIf(not scipy, "scipy required")
class TestCobraIOmat(TestCase, TestCobraIO):
    def setUp(self):
        self.test_model = mini_model
        self.test_file = join(data_directory, "mini.mat")
        self.read_function = io.load_matlab_model
        self.write_function = io.save_matlab_model

    def extra_comparisons(self, model1, model2):
        # MAT does not store gene names
        None


class TestCobraIOjson(TestCase, TestCobraIO):
    def setUp(self):
        self.test_model = mini_model
        self.test_file = join(data_directory, "mini.json")
        self.read_function = io.load_json_model
        self.write_function = io.save_json_model

    def validate(self, filename):
        with open(filename, "r") as infile:
            loaded = io.json.json.load(infile)
        if jsonschema is None:
            warn("jsonschema not installed")
            return
        jsonschema.validate(loaded, io.json.json_schema)


class TestCobraIOPickle(TestCase, TestCobraIO):
    def setUp(self):
        self.test_model = mini_model
        self.test_file = join(data_directory, "mini.pickle")
        self.load = load
        self.dump = dump

        def read_function(filename):
            with open(filename, "rb") as infile:
                return self.load(infile)

        def write_function(model, filename):
            with open(filename, "wb") as outfile:
                self.dump(model, outfile)

        self.read_function = read_function
        self.write_function = write_function


@skipIf(not cPickle, "cPickle required")
class TestCobraIOcPickle(TestCobraIOPickle):
    def setUp(self):
        TestCobraIOPickle.setUp(self)
        self.load = cPickle.load
        self.dump = cPickle.dump


# make a test suite to run all of the tests
loader = TestLoader()
suite = loader.loadTestsFromModule(sys.modules[__name__])


def test_all():
    TextTestRunner(verbosity=2).run(suite)

if __name__ == "__main__":
    test_all()
