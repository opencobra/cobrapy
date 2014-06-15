try:
    import setuptools
except ImportError:
    import ez_setup
    ez_setup.use_setuptools()
from setuptools import setup, find_packages
from sys import argv

import warnings
with warnings.catch_warnings():
    warnings.simplefilter("ignore")
    from cobra.version import get_version

__version = get_version(pep440=True)
setup_kwargs = {}

# for running parallel tests due to a bug in python 2.7.3
# http://bugs.python.org/issue15881#msg170215
try:
    import multiprocessing
except:
    None

# cython is optional for building. The c file can be used directly. However,
# to run sdist, the c file must be generated, which requires cython.
try:
    from Cython.Build import cythonize
except ImportError:
    cythonize = None
    if "sdist" in argv:
        raise Exception("cython required for sdist")

# for building the cglpk solver
try:
    from distutils.extension import Extension
    from distutils.command.build_ext import build_ext
    from os.path import isfile, abspath, dirname
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
    if system() == "Darwin":  # otherwise Mac Clang gives errors
        build_args["extra_compile_args"] = ["-Qunused-arguments"]
    build_args["libraries"] = ["glpk"]
    setup_kwargs["cmdclass"] = {"build_ext": FailBuild}
    # To statically link libglpk to the built extension, add the glpk.h header
    # and static library libglpk.a to the build directory. A static libglpk.a
    # can be built by running configure with the export CLFAGS="-fPIC" and
    # copying the file from src/.libs
    if isfile("libglpk.a"):
        build_args["library_dirs"] = [dirname(abspath("libglpk.a"))]
    if isfile("glpk.h"):
        build_args["include_dirs"] = [dirname(abspath("glpk.h"))]
    # use cython if present, otherwise use c file
    if cythonize:
        ext_modules = cythonize([Extension("cobra.solvers.cglpk",
                ["cobra/solvers/cglpk.pyx"], **build_args)])
    else:
        ext_modules = [Extension("cobra.solvers.cglpk",
                ["cobra/solvers/cglpk.c"], **build_args)]
except:
    ext_modules = None

extras = {
    'parallel': ['pp>=1.6.0'],
    'matlab': ["mlabwrap>=1.1"],
    'sbml': ["python-libsbml-experimental"],
    'array': ["numpy>=1.6", "scipy>=11.0"],
    'display': ["matplotlib", "brewer2mpl", "pandas"]
}

all_extras = set()
for extra in extras.values():
    all_extras.update(extra)
extras["all"] = list(all_extras)

setup(
    name = "cobra",
    version = __version,
    packages = find_packages(exclude=['cobra.oven', 'cobra.oven*']),
    #scripts = [''],
    #put in numpy, scipy, libsbml, and pyglpk
    setup_requires = [],
    #install_requires = ['numpy>=1.6', 'scipy>=0.10'],
    #leave blank because it tries to build scipy/numpy on os x when they are
    #installed by the superpack.  And these are not really essential for core functions.
    install_requires = [],
    extras_require = extras,
    ext_modules = ext_modules,

    package_data = {
         '': ['test/data/*',
              'VERSION',
              'mlab/matlab_scripts/*m']},

    author = "Daniel Robert Hyduke <danielhyduke@gmail.com>, Ali Ebrahim <aebrahim@ucsd.edu>",
    author_email = "danielhyduke@gmail.com",
    description = "COBRApy is a package for constraints-based modeling of biological networks",
    license = "GPL V3.0",
    keywords = "metabolism biology linear programming optimization flux balance analysis fba",
    url = "https://opencobra.github.io/cobrapy",
    test_suite = "cobra.test.suite",
    long_description = "COnstraint-Based Reconstruction and Analysis (COBRA) methods are widely used for genome-scale modeling of metabolic networks in both prokaryotes and eukaryotes.  COBRApy is a constraint-based modeling package that is designed to accomodate the biological complexity of the next generation of COBRA models and provides access to commonly used COBRA methods, such as flux balance analysis, flux variability analysis, and gene deletion analyses.  Through the mlabwrap module it is possible to use COBRApy to call many additional COBRA methods present in the COBRA Toolbox for MATLAB.",
    download_url = 'https://pypi.python.org/pypi/cobra',
    classifiers = ['Development Status :: 5 - Production/Stable',
                   'Environment :: Console',
                   'Intended Audience :: Science/Research',
                   'License :: OSI Approved :: GNU General Public License v3 or later (GPLv3+)',
                   'Operating System :: MacOS :: MacOS X',
                   'Operating System :: Microsoft :: Windows :: Windows 7',
                   'Operating System :: Microsoft :: Windows :: Windows Vista',
                   'Operating System :: Microsoft :: Windows :: Windows XP',
                   'Operating System :: POSIX',
                   'Operating System :: POSIX :: Linux',
                   'Programming Language :: Python :: 2.7',
                   'Programming Language :: Python :: 3.3',
                   'Programming Language :: Python :: 3.4',
                   'Programming Language :: Python :: Implementation :: CPython',
                   'Programming Language :: Python :: Implementation :: Jython',
                   'Topic :: Scientific/Engineering',
                   'Topic :: Scientific/Engineering :: Bio-Informatics'
                   ],
    platforms = "Python >= 2.6 on GNU/Linux, Mac OS X >= 10.7, Microsoft Windows >= 7. \n Jython >= 2.5 on Java >= 1.6",
    **setup_kwargs
    )
