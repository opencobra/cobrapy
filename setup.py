# -*- coding: utf-8 -*-

from os.path import abspath, dirname, isfile, join
from sys import argv
from warnings import warn

from setuptools import setup, find_packages

# cython is optional for building. The c file can be used directly. However,
# for certain functions, the c file must be generated, which requires cython.
try:
    from Cython.Build import cythonize
    from distutils.version import StrictVersion
    import Cython
    try:
        cython_version = StrictVersion(Cython.__version__)
    except ValueError:
        raise ImportError("Cython version not parseable")
    else:
        if cython_version < StrictVersion("0.21"):
            raise ImportError("Cython version too old to use")
except ImportError:
    cythonize = None
    for k in ["sdist", "develop"]:
        if k in argv:
            raise Exception("Cython >= 0.21 required for " + k)

# Begin constructing arguments for building
setup_kwargs = {}

# for building the cglpk solver
try:
    from distutils.extension import Extension
    from distutils.command.build_ext import build_ext
    from os import name
    from platform import system

    class FailBuild(build_ext):
        """allow building of the C extension to fail"""
        def run(self):
            try:
                build_ext.run(self)
            except Exception as e:
                warn(e)

        def build_extension(self, ext):
            try:
                build_ext.build_extension(self, ext)
            except:
                None

    build_args = {}
    setup_kwargs["cmdclass"] = {"build_ext": FailBuild}
    # MAC OS X needs some additional configuration tweaks
    # Build should be run with the python.org python
    # Cython will output C which could generate warnings in clang
    # due to the addition of additional unneeded functions. Because
    # this is a known phenomenon, these warnings are silenced to
    # make other potential warnings which do signal errors stand
    # out.
    if system() == "Darwin":
        build_args["extra_compile_args"] = ["-Wno-unused-function"]

    build_args["libraries"] = ["glpk"]
    # It is possible to statically link libglpk to the built extension. This
    # allows for simplified installation without the need to install libglpk to
    # the system, and is also usueful when installing a particular version of
    # glpk which conflicts with thesystem version. A static libglpk.a can be
    # built by running configure with the export CLFAGS="-fPIC" and copying the
    # file from src/.libs to either the default lib directory or to the build
    # directory. For an example script, see
    # https://gist.github.com/aebrahim/94a2b231d86821f7f225
    include_dirs = []
    library_dirs = []
    if isfile("libglpk.a"):
        library_dirs.append(abspath("."))
    if isfile("glpk.h"):
        include_dirs.append(abspath("."))
    # if the glpk files are not in the current directory attempt to
    # auto-detect their location by finding the location of the glpsol
    # command
    if name == "posix" and len(include_dirs) == 0 and len(library_dirs) == 0:
        from subprocess import check_output
        try:
            glpksol_path = check_output(["which", "glpsol"],
                                        universal_newlines=True).strip()
            glpk_path = abspath(join(dirname(glpksol_path), ".."))
            include_dirs.append(join(glpk_path, "include"))
            library_dirs.append(join(glpk_path, "lib"))
        except Exception as e:
            print('Could not autodetect include and library dirs: ' + str(e))
    if len(include_dirs) > 0:
        build_args["include_dirs"] = include_dirs
    if len(library_dirs) > 0:
        build_args["library_dirs"] = library_dirs
    # use cython if present, otherwise use c file
    if cythonize:
        ext_modules = cythonize([Extension("cobra.solvers.cglpk",
                                           ["cobra/solvers/cglpk.pyx"],
                                           **build_args)],
                                force=True)
    else:
        ext_modules = [Extension("cobra.solvers.cglpk",
                                 ["cobra/solvers/cglpk.c"], **build_args)]
except Exception as e:
    print('Could not build CGLPK: {}'.format(e))
    ext_modules = None

setup_requirements = []
# prevent pytest-runner from being installed on every invocation
if {'pytest', 'test', 'ptr'}.intersection(argv):
    setup_requirements.append("pytest-runner")

extras = {
    'matlab': ["pymatbridge"],
    'sbml': ["python-libsbml", "lxml"],
    'array': ["scipy>=0.11.0"],
    'display': ["matplotlib", "palettable"]
}

all_extras = {'Cython>=0.21'}
for extra in extras.values():
    all_extras.update(extra)
extras["all"] = sorted(list(all_extras))

# If using bdist_wininst, the installer will not get dependencies like
# a setuptools installation does. Therefore, for the one external dependency,
# which is six.py, we can just download it here and include it in the
# installer.

# The file six.py will need to be manually downloaded and placed in the
# same directory as setup.py.
if "bdist_wininst" in argv:
    setup_kwargs["py_modules"] = ["six"]

try:
    with open('README.rst') as handle:
        readme = handle.read()
    with open('INSTALL.rst') as handle:
        install = handle.read()
    setup_kwargs["long_description"] = readme + "\n\n" + install
except:
    setup_kwargs["long_description"] = ''

setup(
    name="cobra",
    version="0.6.0a7",
    packages=find_packages(),
    setup_requires=setup_requirements,
    install_requires=["future", "swiglpk", "optlang>=1.1.4", "ruamel.yaml"
                      "pandas>=0.17.0", "numpy>=1.6", "tabulate"],
    tests_require=["jsonschema > 2.5", "pytest", "pytest-benchmark"],
    extras_require=extras,
    ext_modules=ext_modules,

    package_data={
         '': ['test/data/*',
              'mlab/matlab_scripts/*m']},

    author="Daniel Robert Hyduke <danielhyduke@gmail.com>, "
    "Ali Ebrahim <aebrahim@ucsd.edu>",
    author_email="aebrahim@ucsd.edu",
    description="COBRApy is a package for constraints-based modeling of "
    "biological networks",
    license="LGPL/GPL v2+",
    keywords="metabolism biology linear programming optimization flux"
    " balance analysis fba",
    url="https://opencobra.github.io/cobrapy",
    test_suite="cobra.test.suite",
    download_url='https://pypi.python.org/pypi/cobra',
    classifiers=[
        'Development Status :: 5 - Production/Stable',
        'Environment :: Console',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: GNU Lesser General Public License v2'
            ' or later (LGPLv2+)',
        'License :: OSI Approved :: GNU General Public License v2'
            ' or later (GPLv2+)',
        'Operating System :: OS Independent',
        'Programming Language :: Python :: 2.7',
        'Programming Language :: Python :: 3.4',
        'Programming Language :: Python :: 3.5',
        'Programming Language :: Python :: 3.6',
        'Programming Language :: Cython',
        'Programming Language :: Python :: Implementation :: CPython',
        'Topic :: Scientific/Engineering',
        'Topic :: Scientific/Engineering :: Bio-Informatics'
    ],
    platforms="GNU/Linux, Mac OS X >= 10.7, Microsoft Windows >= 7",
    **setup_kwargs)
