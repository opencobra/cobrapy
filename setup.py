import ez_setup
ez_setup.use_setuptools()
from setuptools import setup, find_packages
from cobra.version import get_version

__version = get_version()
setup_kwargs = {}

# for running parallel tests due to a bug in python 2.7.3
# http://bugs.python.org/issue15881#msg170215
try:
    import multiprocessing
except:
    None

try:
    from Cython.Build import cythonize
    from distutils.extension import Extension
    from distutils.command.build_ext import build_ext
    from os.path import isfile, abspath, dirname
    from platform import system

    class FailBuild(build_ext):
        """allow building of the C extension to fail"""
        def run(self):
            try:
                build_ext.run(self)
            except:
                None

        def build_extension(self, ext):
            try:
                build_ext.build_extension(self, ext)
            except:
                None

    sources = ["cobra/solvers/cglpk.pyx"]
    build_args = {}
    if system() == "Darwin":
        build_args["extra_compile_args"] = ["-Qunused-arguments"]
    build_args["libraries"] = ["glpk"]
    setup_kwargs["cmdclass"] = {"build_ext": FailBuild}
    if isfile("libglpk.a"):
        build_args["library_dirs"] = [dirname(abspath("libglpk.a"))]
    if isfile("glpk.h"):
        build_args["include_dirs"] = [dirname(abspath("glpk.h"))]
    ext_modules = cythonize([Extension("cobra.solvers.cglpk", sources,
                                       **build_args)])
except:
    ext_modules = None

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
    extras_require = {
        'parallel': ['pp>=1.6.0'],
        'matlab': ["mlabwrap>=1.1"],
        'R': ["rpy2>=2.2.2"]
        },
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
    url = "https://github.com/opencobra/cobrapy",
    test_suite = "cobra.test.suite",
    long_description = "COnstraint-Based Reconstruction and Analysis (COBRA) methods are widely used for genome-scale modeling of metabolic networks in both prokaryotes and eukaryotes.  COBRApy is a constraint-based modeling package that is designed to accomodate the biological complexity of the next generation of COBRA models and provides access to commonly used COBRA methods, such as flux balance analysis, flux variability analysis, and gene deletion analyses.  Through the mlabwrap module it is possible to use COBRApy to call many additional COBRA methods present in the COBRA Toolbox for MATLAB.",
    download_url = 'http://sourceforge.net/projects/opencobra/files/python/cobra/' + __version,
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
                   'Programming Language :: Python :: 2.5',
                   'Programming Language :: Python :: 2.6',
                   'Programming Language :: Python :: 2.7',
                   'Programming Language :: Python :: Implementation :: CPython',
                   'Programming Language :: Python :: Implementation :: Jython',
                   'Topic :: Scientific/Engineering',
                   'Topic :: Scientific/Engineering :: Bio-Informatics'
                   ],
    platforms = "Python >= 2.6 on GNU/Linux, Mac OS X >= 10.7, Microsoft Windows >= 7. \n Jython >= 2.5 on Java >= 1.6",
    **setup_kwargs
    )
    
