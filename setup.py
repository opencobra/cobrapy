import ez_setup
ez_setup.use_setuptools()
from setuptools import setup, find_packages
__version = '0.3.0-dev'

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

    package_data = {
         '': ['test/data/*',
              'examples/*py',
              'mlab/matlab_scripts/*m']},

    author = "Daniel Robert Hyduke",
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
    platforms = "Python >= 2.6 on GNU/Linux, Mac OS X >= 10.7, Microsoft Windows >= 7. \n Jython >= 2.5 on Java >= 1.6"

    )
    
