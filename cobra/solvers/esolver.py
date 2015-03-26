from subprocess import check_output, check_call, CalledProcessError
from os import unlink, devnull
from os.path import isfile
from tempfile import NamedTemporaryFile
from fractions import Fraction
from six.moves import zip

from . import cglpk
from .wrappers import *

# detect paths to system calls for esolver and gzip
with open(devnull, "w") as DEVNULL:
    try:
        ESOLVER_COMMAND = check_output(["which", "esolver"],
                                       stderr=DEVNULL).strip()
        __esolver_version__ = check_output(["esolver", "-v"], stderr=DEVNULL)
    except CalledProcessError:
        raise RuntimeError("esolver command not found")
    try:
        GZIP_COMMAND = check_output(["which", "gzip"], stderr=DEVNULL).strip()
    except CalledProcessError:
        raise RuntimeError("gzip command not found")
del DEVNULL

solver_name = "esolver"


class Esolver(cglpk.GLP):
    """contain an LP which will be solved through the QSopt_ex

    The LP is stored using a GLPK object, and written out to an
    LP file which is then solved by the esolver command."""

    def __init__(self, cobra_model=None):
        cglpk.GLP.__init__(self, cobra_model)
        self.solution_filepath = None
        self.basis_filepath = None
        self.rational_solution = False
        self.verbose = False
        self.clean_up = True  # clean up files

    def _clean(self, filename):
        """remove old files"""
        if self.clean_up and filename is not None and isfile(filename):
            unlink(filename)

    def set_parameter(self, parameter_name, value):
        if parameter_name == "GLP":
            raise Exception("can not be set this way")
        if parameter_name == "objective_sense":
            self.set_objective_sense(value)
        if not hasattr(self, parameter_name):
            raise ValueError("Unkonwn parameter '%s'" % parameter_name)
        setattr(self, parameter_name, value)

    def solve_problem(self, **solver_parameters):
        if "objective_sense" in solver_parameters:
            self.set_objective_sense(solver_parameters.pop("objective_sense"))
        for key, value in solver_parameters.items():
            self.set_parameter(key, value)
        # remove the old solution file
        self._clean(self.solution_filepath)
        with NamedTemporaryFile(suffix=".lp", delete=False) as f:
            lp_filepath = f.name
        self.write(lp_filepath)
        existing_basis = self.basis_filepath
        with NamedTemporaryFile(suffix=".bas", delete=False) as f:
            self.basis_filepath = f.name
        with NamedTemporaryFile(suffix=".sol") as f:
            self.solution_filepath = f.name
        command = [ESOLVER_COMMAND, "-b", self.basis_filepath,
                   "-O", self.solution_filepath[:-4]]
        if existing_basis is not None and isfile(existing_basis):
            command.extend(["-B", existing_basis])
        command.extend(["-L", lp_filepath])
        command_kwargs = {}
        if self.verbose:
            print(" ".join(command))
            DEVNULL = None
        else:
            DEVNULL = open(devnull, 'wb')
            command_kwargs["stdout"] = DEVNULL
            command_kwargs["stderr"] = DEVNULL
        try:
            check_call(command, **command_kwargs)
            failed = False
        except CalledProcessError as e:
            failed = True
        if failed:
            self.basis_filepath = existing_basis
            existing_basis = None
            # Sometimes on failure a solution isn't written out
            if not isfile(self.solution_filepath):
                with open(self.solution_filepath, "w") as outfile:
                    outfile.write("=infeasible\n")
        elif isfile(self.solution_filepath + ".gz"):
            # the solution may be written out compressed
            check_call([GZIP_COMMAND, "-d", self.solution_filepath + ".gz"])
        if DEVNULL is not None:
            DEVNULL.close()
        self._clean(lp_filepath)
        self._clean(existing_basis)  # replaced with the new basis

    def get_status(self):
        with open(self.solution_filepath) as infile:
            return infile.readline().split("=")[1].strip().lower()

    def _format(self, value):
        """convert a string value into either a fraction or float"""
        value = Fraction(value)
        return value if self.rational_solution else float(value)

    def get_objective_value(self):
        with open(self.solution_filepath) as infile:
            status = infile.readline().split("=")[1].strip().lower()
            if status != "optimal":
                raise RuntimeError("status not optimal")
            infile.readline()
            return self._format(infile.readline().split("=")[1].strip())

    def format_solution(self, cobra_model):
        m = cobra_model
        solution = m.solution.__class__(None)
        with open(self.solution_filepath) as infile:
            solution.status = infile.readline().split("=")[1].strip().lower()
            if solution.status != "optimal":
                return solution
            infile.readline()
            solution.f = self._format(Fraction(infile.readline()
                                               .split("=")[1].strip()))
            infile.readline()
            value_dict = {}
            for line in infile:
                if line.endswith(":\n"):
                    break
                varname, value = line.split("=")
                value_dict[varname.strip()] = self._format(value.strip())
            dual_dict = {}
            for line in infile:
                if line.endswith(":\n"):
                    break
                varname, value = line.split("=")
                dual_dict[varname.strip()] = self._format(value.strip())
        solution.x = [value_dict.get("x_%d" % (i + 1), 0)
                      for i in range(len(m.reactions))]
        solution.x_dict = {r.id: v for r, v in zip(m.reactions, solution.x)}
        solution.y = [dual_dict.get("r_%d" % (i + 1), 0)
                      for i in range(len(m.metabolites))]
        solution.y_dict = {m.id: v for m, v in zip(m.metabolites, solution.y)}
        return solution


# wrappers for the classmethods at the module level
create_problem = Esolver.create_problem
solve = Esolver.solve
