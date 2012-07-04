import ez_setup
ez_setup.use_setuptools()
from setuptools import setup, find_packages
from pdb import set_trace
setup(
    name = "cobra",
    version = '0.1.0',
    packages = find_packages(exclude=['cobra.internal', 'cobra.oven', 'cobra.db_tools',
                                      'cobra.omics', 'cobra.general',]),
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

         '': ['*.txt', '*.html','LICENSE','README','test/data/*',
              'documentation(/*)+[(html)|(txt)|(png)|(css)|(gif)|(js)]',
              'examples/*py']},

    author = "Daniel Robert Hyduke",
    author_email = "danielhyduke@gmail.com",
    description = "COBRApy is a package for constraints-based modeling of biological networks",
    license = "GPL V3.0",
    keywords = "metabolism biology linear programming optimization flux balance analysis fba",
    url = "http://opencobra.sourceforge.net"
    )
    
